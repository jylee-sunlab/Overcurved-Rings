function res = overcurved_energy_scan(p, varargin)
%OVERCURVED_ENERGY_SCAN
%
%   RES = OVERCURVED_ENERGY_SCAN(P) minimizes the elastic energy of an
%   overcurved ring over a sweep of prescribed overcurving ratios O_p
%   and returns the equilibrium overcurving ratio O_geom
%
%   INPUT
%   -----
%   P is a struct with the following required fields (SI units)
%
%       P.m     lobe-pair count.  Must be 2.
%       P.E     Young's modulus                                      [Pa]
%       P.nu    Poisson's ratio                                      [-]
%       P.I1    second moment of the strong bending channel          [m^4]
%               (the channel that carries the preset curvature)
%       P.I2    second moment of the weak bending channel            [m^4]
%       P.J     Saint-Venant torsional constant                      [m^4]
%       P.A     cross-sectional area                                 [m^2]
%       P.kp    preset (natural) curvature, kappa_p = 1/R_0          [1/m]
%
%   Optional fields
%
%       P.Ns0            (201)   spatial grid points per lobe.  MUST BE ODD.
%       P.nDivTheta0     (100)   scan intervals; the scan has nDivTheta0+1
%                                points spanning O_p in [1, 2*m-1].
%       P.absReg         (1e-10) regularization of |sin(phi)cos(t cos phi)|
%       P.tanMargin      (1e-8)  margin keeping t*cos(phi) away from pi/2
%       P.tauSign        (+1)    sign convention of the geometric torsion
%       P.qBound         (6.0)   bound on the map generator q
%       P.fdRelStepOmega (1e-6)  relative FD step for d/d(omega)
%       P.fdRelStepQ     (1e-6)  relative FD step for d/dq
%       P.nProbeEvery    (4)     probe the opposite basin every N-th step
%       P.dwellAfterB    (6)     steps to keep dual mode after the probe wins
%
%   NAME-VALUE OPTIONS
%   ------------------
%       'Verbose'   (true)  print per-step progress
%       'PrintEvery'(1)     print every N-th step when Verbose is true
%
%   OUTPUT
%   ------
%   RES is a struct with fields
%
%       RES.table    table, one row per scan point, with columns
%                    Op, theta0, L0, omega, phi, ell, Oeq,
%                    Ub, Ut, Ua, Utotal, exitflag
%       RES.fields   struct of cell arrays, one cell per scan point:
%                    s0, theta_m, t, dt_ds0, lambda, q
%       RES.diag     struct of per-step solver diagnostics
%       RES.p        the validated and defaulted input struct, with the
%                    derived quantities G, EI1, EI2, GJ, EA appended
%       RES.meta     code version, MATLAB release, timestamp, wall time
%
%   REQUIREMENTS
%   ------------
%   MATLAB with the Optimization Toolbox (fmincon).  No other toolboxes.

CODE_VERSION = '1.0.0';

opts = parse_options(varargin{:});
p    = validate_and_default(p);

%% grid
theta0_min  = pi / p.m;
theta0_max  = 2*pi - pi / p.m;
theta0_list = linspace(theta0_min, theta0_max, p.nDivTheta0 + 1).';
nPts        = numel(theta0_list);
Op_list     = p.m * theta0_list / pi;

%% storage
nan_col = nan(nPts, 1);

L0_c     = nan_col;   omega_c = nan_col;   phi_c    = nan_col;
ell_c    = nan_col;   Oeq_c   = nan_col;   Ub_c     = nan_col;
Ut_c     = nan_col;   Ua_c    = nan_col;   Utot_c   = nan_col;
exitfl_c = nan_col;

fields = struct( ...
    's0',      {cell(nPts,1)}, ...
    'theta_m', {cell(nPts,1)}, ...
    't',       {cell(nPts,1)}, ...
    'dt_ds0',  {cell(nPts,1)}, ...
    'lambda',  {cell(nPts,1)}, ...
    'q',       {cell(nPts,1)});

diag = struct( ...
    'startWinner',   nan_col, ...
    'fvalA',         nan_col, ...
    'fvalB',         nan_col, ...
    'exitA',         nan_col, ...
    'exitB',         nan_col, ...
    'ranB',          false(nPts,1), ...
    'iterations',    nan_col, ...
    'funcCount',     nan_col, ...
    'constrViolation', nan_col, ...
    'firstOrderOpt', nan_col);

fminOpt = make_fmincon_options();

prevSolA       = [];
prevSolB       = [];
dwellRemaining = 0;

if opts.Verbose
    fprintf('overcurved_energy_scan %s\n', CODE_VERSION);
    fprintf('  Op in [%.4f, %.4f], %d points, Ns0 = %d\n', ...
            Op_list(1), Op_list(end), nPts, p.Ns0);
    fprintf('  EI1 = %.4e, EI2 = %.4e, GJ = %.4e, EA = %.4e\n', ...
            p.EI1, p.EI2, p.GJ, p.EA);
end

tStart = tic;

%% scan
for i = 1:nPts

    theta0 = theta0_list(i);
    L0     = theta0 / p.kp;

    s0  = linspace(0, L0, p.Ns0).';
    ds0 = s0(2) - s0(1);
    ic  = (p.Ns0 + 1) / 2;

    cache = make_qmap_cache(s0, ds0);

    idxThetaFree = 2:(p.Ns0-1);
    nThetaFree   = numel(idxThetaFree);
    nQLeft       = ic;

    lb = [p.omega_min + 1e-8;
          log(1e-8);
          -inf(nThetaFree,1);
          -p.qBound * ones(nQLeft,1)];

    ub = [p.omega_max - 1e-8;
          log(1e+8);
           inf(nThetaFree,1);
           p.qBound * ones(nQLeft,1)];

    objF = @(x) objective_grad(x, s0, idxThetaFree, p, cache, lb, ub);
    nlcF = @(x) nonlinear_constraints(x, s0, ds0, idxThetaFree, p);

    % ---- Start A: warm continuation from the previous scan point -----
    if isempty(prevSolA)
        omega0_A = clampd(theta0, p.omega_min + 1e-6, p.omega_max - 1e-6);
        ell0_A   = 1 / p.kp;
        theta_gA = zeros(p.Ns0, 1);
        qLeft_gA = zeros(nQLeft, 1);
    else
        omega0_A = clampd(prevSolA.omega, p.omega_min + 1e-6, p.omega_max - 1e-6);
        ell0_A   = prevSolA.ell;
        theta_gA = prevSolA.theta;
        qLeft_gA = prevSolA.qfull(1:ic);
    end

    x0_A = pack_design(omega0_A, ell0_A, theta_gA, qLeft_gA, idxThetaFree);

    try
        [xA, fvalA, exitA, outputA] = fmincon( ...
            objF, x0_A, [], [], [], [], lb, ub, nlcF, fminOpt);
    catch ME
        warning('overcurved:startAFailed', ...
                'Step %d: Start A failed (%s).', i, ME.message);
        xA = x0_A;  fvalA = inf;  exitA = -99;  outputA = struct();
    end

    % ---- Decide whether to probe the opposite basin ------------------
    runB = (i == 1) || (mod(i, p.nProbeEvery) == 0) || (dwellRemaining > 0);

    if runB
        if isempty(prevSolA)
            omega_warm = theta0;
        else
            omega_warm = prevSolA.omega;
        end

        if omega_warm < pi
            omega_target_B = 3*pi/2;
            ell_target_B   = 1 / (3 * p.kp);
            target_high    = true;
        else
            omega_target_B = pi/2;
            ell_target_B   = 1 / p.kp;
            target_high    = false;
        end

        canWarmB = false;
        if ~isempty(prevSolB)
            if ( target_high && prevSolB.omega >  pi) || ...
               (~target_high && prevSolB.omega <= pi)
                canWarmB = true;
            end
        end

        if canWarmB
            omega0_B = clampd(prevSolB.omega, p.omega_min + 1e-6, p.omega_max - 1e-6);
            ell0_B   = prevSolB.ell;
            theta_gB = prevSolB.theta;
            qLeft_gB = prevSolB.qfull(1:ic);
        else
            omega0_B = clampd(omega_target_B, p.omega_min + 1e-6, p.omega_max - 1e-6);
            ell0_B   = ell_target_B;
            theta_gB = zeros(p.Ns0, 1);
            qLeft_gB = zeros(nQLeft, 1);
        end

        x0_B = pack_design(omega0_B, ell0_B, theta_gB, qLeft_gB, idxThetaFree);

        try
            [xB, fvalB, exitB, outputB] = fmincon( ...
                objF, x0_B, [], [], [], [], lb, ub, nlcF, fminOpt);
        catch ME
            warning('overcurved:startBFailed', ...
                    'Step %d: Start B failed (%s).', i, ME.message);
            xB = x0_B;  fvalB = inf;  exitB = -99;  outputB = struct();
        end
    else
        xB = [];  fvalB = NaN;  exitB = NaN;  outputB = struct();
    end

    % ---- Pick the lower of the two ----------------------------------
    validA = isfinite(fvalA) && fvalA < 1e29;
    validB = runB && isfinite(fvalB) && fvalB < 1e29;

    if validA && validB
        if fvalB < fvalA - 1e-12
            xopt = xB;  exitflag = exitB;  output = outputB;  winner = 2;
        else
            xopt = xA;  exitflag = exitA;  output = outputA;  winner = 1;
        end
    elseif validA
        xopt = xA;  exitflag = exitA;  output = outputA;  winner = 1;
    elseif validB
        xopt = xB;  exitflag = exitB;  output = outputB;  winner = 2;
    else
        xopt = xA;  exitflag = exitA;  output = outputA;  winner = 1;
    end

    if winner == 2
        dwellRemaining = p.dwellAfterB;
    else
        dwellRemaining = max(0, dwellRemaining - 1);
    end

    % ---- Reconstruct the winning state ------------------------------
    [omega, ell, theta, qfull] = unpack_design(xopt, p.Ns0, idxThetaFree);
    [t, dt_ds0, ~]             = build_t_from_qmap(omega, qfull, s0);

    phi          = phi_from_omega(omega, p.m);
    [Utot, comp] = energy_core(theta, t, dt_ds0, omega, s0, p, ell, cache);

    Oeq = (2*p.m/pi) * tan(phi) * sin(0.5*omega*cos(phi));

    % ---- Store -------------------------------------------------------
    L0_c(i)     = L0;
    omega_c(i)  = omega;
    phi_c(i)    = phi;
    ell_c(i)    = ell;
    Oeq_c(i)    = Oeq;
    Ub_c(i)     = comp.Ub_total;
    Ut_c(i)     = comp.Ut_total;
    Ua_c(i)     = comp.Ua_total;
    Utot_c(i)   = Utot;
    exitfl_c(i) = exitflag;

    fields.s0{i}      = s0;
    fields.theta_m{i} = theta;
    fields.t{i}       = t;
    fields.dt_ds0{i}  = dt_ds0;
    fields.lambda{i}  = comp.lambda;
    fields.q{i}       = qfull;

    diag.startWinner(i)     = winner;
    diag.fvalA(i)           = fvalA;
    diag.fvalB(i)           = fvalB;
    diag.exitA(i)           = exitA;
    diag.exitB(i)           = exitB;
    diag.ranB(i)            = runB;
    diag.iterations(i)      = get_field(output, 'iterations',      NaN);
    diag.funcCount(i)       = get_field(output, 'funcCount',       NaN);
    diag.constrViolation(i) = get_field(output, 'constrviolation', NaN);
    diag.firstOrderOpt(i)   = get_field(output, 'firstorderopt',   NaN);

    prevSolA.omega = omega;
    prevSolA.ell   = ell;
    prevSolA.theta = theta;
    prevSolA.qfull = qfull;

    if runB && validB
        [omB, ellB, thB, qfB] = unpack_design(xB, p.Ns0, idxThetaFree);
        prevSolB.omega = omB;
        prevSolB.ell   = ellB;
        prevSolB.theta = thB;
        prevSolB.qfull = qfB;
    end

    % ---- Progress ----------------------------------------------------
    if opts.Verbose && (mod(i, opts.PrintEvery) == 0 || i == 1 || i == nPts)
        if winner == 1
            who = 'A';
        else
            who = 'B';
        end
        fprintf('  [%4d/%d] Op=%6.4f  Oeq=%7.5f  U=%.6e  win=%s  (%.1fs)\n', ...
                i, nPts, Op_list(i), Oeq, Utot, who, toc(tStart));
    end
end

elapsed = toc(tStart);

%% output
res.table = table(Op_list, theta0_list, L0_c, omega_c, phi_c, ell_c, ...
                  Oeq_c, Ub_c, Ut_c, Ua_c, Utot_c, exitfl_c, ...
    'VariableNames', {'Op','theta0','L0','omega','phi','ell', ...
                      'Oeq','Ub','Ut','Ua','Utotal','exitflag'});

res.fields = fields;
res.diag   = diag;
res.p      = p;

res.meta = struct( ...
    'codeVersion',    CODE_VERSION, ...
    'matlabRelease',  version('-release'), ...
    'timestamp',      char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'elapsedSeconds', elapsed);

if opts.Verbose
    fprintf('Done in %.1f s.  Oeq departs from 1 at Op = %s.\n', ...
            elapsed, departure_string(Op_list, Oeq_c));
end
end


%% ========================================================================
%  Input handling
%% ========================================================================
function opts = parse_options(varargin)
opts.Verbose    = true;
opts.PrintEvery = 1;

if mod(numel(varargin), 2) ~= 0
    error('overcurved:badOptions', ...
          'Name-value options must come in pairs.');
end

for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~(ischar(name) || isstring(name))
        error('overcurved:badOptions', 'Option names must be text.');
    end
    switch lower(char(name))
        case 'verbose'
            opts.Verbose = logical(varargin{k+1});
        case 'printevery'
            opts.PrintEvery = max(1, round(varargin{k+1}));
        otherwise
            error('overcurved:badOptions', ...
                  'Unknown option ''%s''.', char(name));
    end
end
end


function p = validate_and_default(p)
if ~isstruct(p) || ~isscalar(p)
    error('overcurved:badInput', 'P must be a scalar struct.');
end

required = {'m','E','nu','I1','I2','J','A','kp'};
missing  = required(~isfield(p, required));
if ~isempty(missing)
    error('overcurved:missingField', ...
          'Missing required field(s): %s.', strjoin(missing, ', '));
end

for k = 1:numel(required)
    v = p.(required{k});
    if ~(isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v))
        error('overcurved:badInput', ...
              'P.%s must be a finite real scalar.', required{k});
    end
end

% --- scope ----------------------------------------------------------
if p.m ~= 2
    error('overcurved:unsupportedM', ...
         ['Only m = 2 is supported (got m = %g).  For m >= 3 the ' ...
          'torsion field has a discontinuous first derivative at the ' ...
          'lobe junctions, and the basin targets used by the solver ' ...
          'are specific to m = 2.'], p.m);
end

% --- positivity -----------------------------------------------------
positive = {'E','I1','I2','J','A','kp'};
for k = 1:numel(positive)
    if p.(positive{k}) <= 0
        error('overcurved:badInput', ...
              'P.%s must be positive (got %g).', ...
              positive{k}, p.(positive{k}));
    end
end

if p.nu <= -1 || p.nu >= 0.5
    error('overcurved:badInput', ...
          'P.nu must lie in (-1, 0.5) (got %g).', p.nu);
end

% --- convention check -----------------------------------------------
if p.I2 > p.I1
    warning('overcurved:weakAxisFirst', ...
           ['P.I2 (%g) exceeds P.I1 (%g).  I1 is the channel carrying ' ...
            'the preset curvature; the paper takes it to be the strong ' ...
            'axis.  Check that the two moments are not swapped.'], ...
            p.I2, p.I1);
end

% --- numerical defaults ---------------------------------------------
p = default_field(p, 'Ns0',            201);
p = default_field(p, 'nDivTheta0',     100);
p = default_field(p, 'absReg',         1e-10);
p = default_field(p, 'tanMargin',      1e-8);
p = default_field(p, 'tauSign',        +1);
p = default_field(p, 'qBound',         6.0);
p = default_field(p, 'fdRelStepOmega', 1e-6);
p = default_field(p, 'fdRelStepQ',     1e-6);
p = default_field(p, 'nProbeEvery',    4);
p = default_field(p, 'dwellAfterB',    6);

if mod(p.Ns0, 2) == 0 || p.Ns0 < 21
    error('overcurved:badGrid', ...
         ['P.Ns0 must be an odd integer >= 21 (got %g).  The solver ' ...
          'imposes lobe symmetry about the midpoint index ' ...
          '(Ns0+1)/2.'], p.Ns0);
end

if p.nDivTheta0 < 1 || mod(p.nDivTheta0, 1) ~= 0
    error('overcurved:badGrid', ...
          'P.nDivTheta0 must be a positive integer (got %g).', ...
          p.nDivTheta0);
end

if p.Ns0 >= 801
    warning('overcurved:expensiveGrid', ...
           ['Ns0 = %d is expensive (order of minutes per scan point) ' ...
            'and will not converge from a cold start.  Keep the full ' ...
            'warm continuation from Op = 1.'], p.Ns0);
end

% --- derived ---------------------------------------------------------
p.G         = p.E / (2*(1 + p.nu));
p.EI1       = p.E * p.I1;
p.EI2       = p.E * p.I2;
p.GJ        = p.G * p.J;
p.EA        = p.E * p.A;
p.omega_min = pi / p.m;
p.omega_max = 2*pi - pi / p.m;
end


function p = default_field(p, name, value)
if ~isfield(p, name) || isempty(p.(name))
    p.(name) = value;
end
end


function opt = make_fmincon_options()
opt = optimoptions('fmincon', ...
    'Algorithm',                  'sqp', ...
    'Display',                    'off', ...
    'MaxIterations',              1500, ...
    'MaxFunctionEvaluations',     1.0e6, ...
    'StepTolerance',              1e-10, ...
    'ConstraintTolerance',        1e-7, ...
    'OptimalityTolerance',        1e-7, ...
    'SpecifyObjectiveGradient',   true, ...
    'SpecifyConstraintGradient',  false);
end


%% ========================================================================
%  Design-vector packing
%% ========================================================================
function x = pack_design(omega, ell, theta, qLeft, idxThetaFree)
x = [omega; log(ell); theta(idxThetaFree); qLeft(:)];
end


function [omega, ell, theta, qfull] = unpack_design(x, Ns0, idxThetaFree)
omega = x(1);
ell   = exp(x(2));

nThetaFree = numel(idxThetaFree);
thetaFree  = x(3 : 2+nThetaFree);
qLeft      = x(3+nThetaFree : end);

theta               = zeros(Ns0, 1);
theta(1)            = 0;
theta(end)          = 0;
theta(idxThetaFree) = thetaFree;

qfull = mirror_q(qLeft, Ns0);
end


function qfull = mirror_q(qLeft, Ns0)
ic = (Ns0 + 1) / 2;

if numel(qLeft) ~= ic
    error('overcurved:qLength', ...
          'qLeft must have (Ns0+1)/2 = %d entries.', ic);
end

qfull = [qLeft(:); qLeft(ic-1:-1:1)];
end


%% ========================================================================
%  Kinematics
%% ========================================================================
function [t, dt_ds0, H] = build_t_from_qmap(omega, qfull, s0)
Ns0 = numel(s0);
ic  = (Ns0 + 1) / 2;
L0  = s0(end) - s0(1);

if L0 <= 0
    error('overcurved:badGrid', 'L0 must be positive.');
end

u     = (s0 - s0(1)) / L0;
uLeft = u(1:ic);

qLeft = qfull(1:ic);
wLeft = exp(qLeft);

WLeft   = cumtrapz(uLeft, wLeft);
halfInt = WLeft(end);

if ~isfinite(halfInt) || halfInt <= 0
    H      = nan(size(s0));
    t      = nan(size(s0));
    dt_ds0 = nan(size(s0));
    return;
end

HLeft    = 0.5 * WLeft / halfInt;
dHduLeft = 0.5 * wLeft / halfInt;

H    = [HLeft;    1 - HLeft(ic-1:-1:1)];
dHdu = [dHduLeft; dHduLeft(ic-1:-1:1)];

H(ic) = 0.5;

t      = omega * (H - 0.5);
dt_ds0 = (omega / L0) * dHdu;
end


function phi = phi_from_omega(omega, m)
val = sin(pi/(2*m)) ./ sin(omega/2);
val = min(max(val, -1), 1);

phi = asin(val);
phi = min(max(phi, 1e-8), pi/2 - 1e-8);
end


function cache = make_qmap_cache(s0, ds0)
n = numel(s0);
e = ones(n, 1);

D = spdiags([-0.5*e, 0.5*e], [-1, 1], n, n) / ds0;

D(1,:) = 0;
D(1,1) = -1/ds0;
D(1,2) =  1/ds0;

D(end,:)     = 0;
D(end,end-1) = -1/ds0;
D(end,end)   =  1/ds0;

wtrap      = ds0 * ones(n, 1);
wtrap(1)   = 0.5 * ds0;
wtrap(end) = 0.5 * ds0;

cache.D     = sparse(D);
cache.wtrap = wtrap;
end


function dfdx = first_derivative(f, dx)
n    = numel(f);
dfdx = zeros(size(f));

if n < 3
    return;
end

dfdx(1)       = (f(2) - f(1)) / dx;
dfdx(end)     = (f(end) - f(end-1)) / dx;
dfdx(2:end-1) = (f(3:end) - f(1:end-2)) / (2*dx);
end


%% ========================================================================
%  Energy
%% ========================================================================
function [Utot, comp] = energy_core(theta, t, dt_ds0, omega, s0, p, ell, cache)
phi = phi_from_omega(omega, p.m);

if ~isfinite(phi) || ~isfinite(ell) || ell <= 0 || phi <= 0 || phi >= pi/2
    Utot = 1e30;
    comp = empty_comp();
    return;
end

theta  = theta(:);
t      = t(:);
dt_ds0 = dt_ds0(:);

if any(~isfinite(dt_ds0)) || any(dt_ds0 <= 0)
    Utot = 1e30;
    comp = empty_comp();
    return;
end

D = cache.D;
w = cache.wtrap(:);

dtheta = D * theta;

g_raw = sin(phi) .* cos(t .* cos(phi));
g_abs = sqrt(g_raw.^2 + p.absReg^2);

h = g_abs .* dt_ds0;

tau_g = (p.tauSign / ell) .* tan(t .* cos(phi));
tau_m = tau_g + dtheta ./ (ell .* h);

kappa_m1 =  (1/ell) .* cos(theta);
kappa_m2 = -(1/ell) .* sin(theta);

lambda = ell .* h;
eps0   = lambda - 1;

bend_density = p.I1 .* (kappa_m1 - p.kp).^2 + p.I2 .* (kappa_m2).^2;
tors_density = tau_m.^2;
axia_density = eps0.^2;

Ub_total = p.m * p.E       * sum(w .* bend_density);
Ut_total = p.m * p.G * p.J * sum(w .* tors_density);
Ua_total = p.m * p.E * p.A * sum(w .* axia_density);

Utot = Ub_total + Ut_total + Ua_total;

if ~isfinite(Utot) || ~isreal(Utot)
    Utot = 1e30;
    comp = empty_comp();
    return;
end

comp.phi      = phi;
comp.dtheta   = dtheta;
comp.g_raw    = g_raw;
comp.g_abs    = g_abs;
comp.h        = h;
comp.tau_g    = tau_g;
comp.tau_m    = tau_m;
comp.kappa_m1 = kappa_m1;
comp.kappa_m2 = kappa_m2;
comp.lambda   = lambda;
comp.eps0     = eps0;

comp.Ub_total = Ub_total;
comp.Ut_total = Ut_total;
comp.Ua_total = Ua_total;
end


function comp = empty_comp()
comp.phi      = nan;
comp.dtheta   = nan;
comp.g_raw    = nan;
comp.g_abs    = nan;
comp.h        = nan;
comp.tau_g    = nan;
comp.tau_m    = nan;
comp.kappa_m1 = nan;
comp.kappa_m2 = nan;
comp.lambda   = nan;
comp.eps0     = nan;

comp.Ub_total = nan;
comp.Ut_total = nan;
comp.Ua_total = nan;
end


%% ========================================================================
%  Objective, gradient, constraints
%% ========================================================================
function [U, grad] = objective_grad(x, s0, idxThetaFree, p, cache, lb, ub)
[omega, ell, theta, qfull] = unpack_design(x, numel(s0), idxThetaFree);
[t, dt_ds0, ~]             = build_t_from_qmap(omega, qfull, s0);

[U, comp] = energy_core(theta, t, dt_ds0, omega, s0, p, ell, cache);

if nargout < 2
    return;
end

grad = zeros(size(x));

if ~isfinite(U) || U >= 1e29
    return;
end

nThetaFree = numel(idxThetaFree);
qStart     = 3 + nThetaFree;

w = cache.wtrap(:);
D = cache.D;

theta = theta(:);
ell2  = ell^2;

k1 = comp.kappa_m1;
k2 = comp.kappa_m2;

% --- d/d theta ------------------------------------------------------
db_dtheta = 2*p.I1 .* (k1 - p.kp) .* (-sin(theta)/ell) ...
          + 2*p.I2 .* (k2)        .* (-cos(theta)/ell);

dUb_dtheta = p.m * p.E * (w .* db_dtheta);

h_theta = ell .* comp.h;
coeff_t = 2 * p.m * p.G * p.J * (w .* comp.tau_m ./ h_theta);

dUt_dtheta = D.' * coeff_t;

g_theta_full = dUb_dtheta + dUt_dtheta;

grad(3 : 2+nThetaFree) = g_theta_full(idxThetaFree);

% --- d/d log(ell) ---------------------------------------------------
db_dell = 2*p.I1 .* (k1 - p.kp) .* (-cos(theta)/ell2) ...
        + 2*p.I2 .* (k2)        .* ( sin(theta)/ell2);

dUb_dell = p.m * p.E * sum(w .* db_dell);
dUt_dell = -2 * comp.Ut_total / ell;
dUa_dell =  2 * p.m * p.E * p.A * sum(w .* comp.eps0 .* comp.h);

dU_dell = dUb_dell + dUt_dell + dUa_dell;

grad(2) = ell * dU_dell;

% --- d/d omega and d/dq (finite differences) ------------------------
valueOnly = @(xx) objective_value(xx, s0, idxThetaFree, p, cache);

grad(1) = fd_component(x, 1, lb, ub, valueOnly, p.fdRelStepOmega);

for k = qStart:numel(x)
    grad(k) = fd_component(x, k, lb, ub, valueOnly, p.fdRelStepQ);
end
end


function U = objective_value(x, s0, idxThetaFree, p, cache)
[omega, ell, theta, qfull] = unpack_design(x, numel(s0), idxThetaFree);
[t, dt_ds0, ~]             = build_t_from_qmap(omega, qfull, s0);

[U, ~] = energy_core(theta, t, dt_ds0, omega, s0, p, ell, cache);
end


function g = fd_component(x, idx, lb, ub, fun, relStep)
xj = x(idx);
h  = relStep * max(1, abs(xj));

if ~isfinite(h) || h <= 0
    h = relStep;
end

hMaxPlus  = ub(idx) - xj;
hMaxMinus = xj - lb(idx);

if (h <= hMaxPlus) && (h <= hMaxMinus)
    xp = x;  xp(idx) = xj + h;
    xm = x;  xm(idx) = xj - h;
    g  = (fun(xp) - fun(xm)) / (2*h);
    return;
end

if h <= hMaxPlus
    xp = x;  xp(idx) = xj + h;
    g  = (fun(xp) - fun(x)) / h;
    return;
end

if h <= hMaxMinus
    xm = x;  xm(idx) = xj - h;
    g  = (fun(x) - fun(xm)) / h;
    return;
end

g = 0;
end


function [c, ceq] = nonlinear_constraints(x, s0, ds0, idxThetaFree, p)
[omega, ~, theta, qfull] = unpack_design(x, numel(s0), idxThetaFree);
[t, ~, ~]                = build_t_from_qmap(omega, qfull, s0);

phi = phi_from_omega(omega, p.m);

dtheta = first_derivative(theta, ds0);
ic     = (numel(s0) + 1) / 2;
ceq    = dtheta(ic);

c = abs(t .* cos(phi)) - (pi/2 - p.tanMargin);
c = c(:);
end


%% ========================================================================
%  Helpers
%% ========================================================================
function y = clampd(x, lo, hi)
y = min(max(x, lo), hi);
end


function val = get_field(s, name, defaultVal)
val = defaultVal;

if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
end
end


function str = departure_string(Op, Oeq)
idx = find(Oeq > 1 + 1e-4, 1, 'first');

if isempty(idx)
    str = 'none in range';
else
    str = sprintf('%.4f', Op(idx));
end
end
