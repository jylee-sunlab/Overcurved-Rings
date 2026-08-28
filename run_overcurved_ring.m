%% Overcurved ring analysis
clear;  clc;


%% ===== 1) USER INPUT =====================================================
p = struct();

% --- Ring symmetry -------------------------------------------------------
p.m  = 2;                % lobe-pair count.  Only m = 2 is supported.

% --- Material ------------------------------------------------------------
p.E  = 2e9;              % Young's modulus [Pa]
p.nu = 0.41;             % Poisson's ratio [1]

% --- Cross-section -------------------------------------------------------
% Supply the four section constants in SI units.
%   I1 : second moment of area (strong axis, I1 >= I2)
%   I2 : second moment of area (weak axis)
%   J  : Saint-Venant torsional constant
%   A  : cross-sectional area
%
% Example: a solid rectangle, 4 mm by 3 mm (Case 6 in the paper):
%   I1 = 3e-3*(4e-3)^3/12,   I2 = 4e-3*(3e-3)^3/12,   A = 4e-3*3e-3,
%   J  = exact Saint-Venant torsional constant.
p.I1 = 1.600000e-11;     % [m^4]
p.I2 = 9.000000e-12;     % [m^4]
p.J  = 1.948939e-11;     % [m^4]
p.A  = 1.200000e-05;     % [m^2]

% --- Preset curvature ----------------------------------------------------
p.kp = 16;               % intrinsic curvature [1/m]

% --- Scan and grid resolution -------------------------------------------
p.Ns0        = 201;      % grid points per lobe.  MUST BE ODD.
p.nDivTheta0 = 200;      % scan intervals. O_p spans [1, 2m-1] in nDivTheta0+1 points

% --- Output --------------------------------------------------------------
outputBaseName = 'rect_4x3_kp16';   % ASCII, no spaces; used for the CSV names
outputFolder   = 'results';         % created if missing; '' writes to pwd
saveFields     = false;             % also write theta_m(s0), t(s0), lambda(s0)
drawShapes     = true;              % also draw the ring shapes of this scan


%% ===== 2) SOLVE ==========================================================
res = overcurved_energy_scan(p);


%% ===== 3) SAVE ===========================================================
save_scan_results(res, outputBaseName, 'Folder', outputFolder, 'SaveFields', saveFields);
matFile = fullfile(outputFolder, [outputBaseName '_results.mat']);
save(matFile, 'res');


%% ===== 4) PLOT ===========================================================
fig = plot_scan_results(res);
plotBase = fullfile(outputFolder, [outputBaseName '_scan']);
savefig(fig, [plotBase '.fig']);
exportgraphics(fig, [plotBase '.png'], 'BackgroundColor', 'white');


%% ===== 5) SHAPES (optional) ==============================================
if drawShapes
    Og = res.table.Oeq(res.table.Oeq > 1 + 1e-4 & res.table.Oeq < 2.996);
    if ~isempty(Og)
        pick = unique(round(linspace(1, numel(Og), min(10, numel(Og)))));
        plot_shapes('Ogeom', Og(pick), ...
                    'OutDir', fullfile(outputFolder, 'shapes'));
    end
end

