function [figs, info] = plot_shapes(varargin)

    opt = parse_args(varargin{:});

    m        = 2;
    O_SWITCH = 1.2;
    O_MAX    = 2*m - 1;

    if isempty(opt.Ogeom)
        Og = linspace(opt.Range(1), opt.Range(2), opt.Count);
    else
        Og = opt.Ogeom(:).';
    end

    if any(Og <= 1) || any(Og >= O_MAX)
        error('overcurved:badOgeom', ...
              ['O_geom must lie strictly between 1 and %g.  ' ...
               'The closed-form shape is singular at both ends.'], O_MAX);
    end

    if opt.Save && ~exist(opt.OutDir, 'dir')
        mkdir(opt.OutDir);
    end

    nO   = numel(Og);
    figs = gobjects(0, 1);
    C    = cell(nO, 1);
    phi  = zeros(nO, 1);  omega = zeros(nO, 1);
    delta = zeros(nO, 1); ysign = zeros(nO, 1);
    gap  = zeros(nO, 1);  arclen = zeros(nO, 1);

    for k = 1:nO
        [C{k}, phi(k), omega(k), delta(k), ysign(k), gap(k), arclen(k)] = ...
            lobe_curves(Og(k), m, O_SWITCH);
    end

    info = table(Og(:), phi, omega, delta, ysign, arclen, gap, ...
        'VariableNames', {'Ogeom','phi','omega','delta','ysign', ...
                          'raw_arclength','closure_gap'});

    fprintf('\n%-9s %-9s %-9s %-9s %-6s %-13s %s\n', ...
            'O_geom', 'phi', 'omega', 'delta', 'ysign', 'raw arclen', ...
            'closure gap');
    for k = 1:nO
        fprintf('%-9.4f %-9.5f %-9.5f %-9.5f %-6s %-13.6f %.2e\n', ...
                Og(k), phi(k), omega(k), delta(k), ...
                sprintf('%+d', ysign(k)), arclen(k), gap(k));
    end
    if any(gap > 1e-6)
        warning('overcurved:shapeClosure', ...
                ['The ring does not close to 1e-6 for %d of %d shapes ' ...
                 '(worst %.2e).  Check the O_geom values.'], ...
                nnz(gap > 1e-6), nO, max(gap));
    end

    for k = 1:nO
        stateFig = figure('Color', 'w', 'Visible', 'off', ...
                          'Position', [100 100 600 600], ...
                          'Name', sprintf('O_geom = %.4f', Og(k)));
        ax = axes(stateFig);
        draw_one(ax, C{k}, opt);

        if opt.Save
            tag = sprintf('shape_Og_%s', ...
                          strrep(sprintf('%.4f', Og(k)), '.', 'p'));
            base = fullfile(opt.OutDir, tag);
            export_png(stateFig, [base '.png'], opt.Resolution);
            savefig(stateFig, [base '.fig']);
            fprintf('  wrote %s.png and .fig\n', tag);
        end

        close(stateFig);
    end

    if opt.Panel
        nc = min(nO, 5);
        nr = ceil(nO/nc);

        panelFig = figure('Color', 'w', ...
                          'Position', [60 60 min(340*nc, 1800) 360*nr], ...
                          'Name', 'Shape progression');
        tl = tiledlayout(panelFig, nr, nc, ...
                         'TileSpacing', 'compact', 'Padding', 'compact');

        for k = 1:nO
            ax = nexttile(tl);
            draw_one(ax, C{k}, opt);
            title(ax, sprintf('$O_{\\mathrm{geom}} = %.3f$', Og(k)), ...
                  'Interpreter', 'latex', 'FontSize', 11, ...
                  'FontWeight', 'normal');
        end

        if opt.Save
            base = fullfile(opt.OutDir, 'shape_progression');
            export_png(panelFig, [base '.png'], opt.Resolution);
            savefig(panelFig, [base '.fig']);
            fprintf('  wrote shape_progression.png and .fig\n');
        end

        figs = panelFig;
    end

    if opt.Save
        fprintf('%d shapes written to %s\n', nO, ...
                fullfile(pwd_if_relative(opt.OutDir)));
    end
end

function draw_one(ax, C, opt)

    axes(ax);
    hold(ax, 'on');
    for j = 1:numel(C)
        lw = opt.LineWidth;
        if opt.HighlightFirstLobe && j == 1
            lw = opt.LineWidth * 1.3;
        end
        plot3(ax, C{j}(1,:), C{j}(2,:), C{j}(3,:), '-', ...
              'Color', [0 0 0], 'LineWidth', lw);
    end
    view(ax, opt.View(1), opt.View(2));
    axis(ax, 'equal');
    axis(ax, 'off');
    grid(ax, 'off');
    box(ax, 'off');
    set(ax, 'Clipping', 'off');
    if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
        ax.Toolbar.Visible = 'off';
    end
    hold(ax, 'off');
end

function [C, phi, omega, delta, ysign, gap, arclen] = ...
                                        lobe_curves(Og, m, O_switch)

    opt = optimoptions('fsolve', 'Display', 'off', ...
                       'TolFun', 1e-12, 'TolX', 1e-12);

    fun = @(x) [ (2*m/pi)*tan(x(1))*sin((x(2)/2)*cos(x(1))) - Og
                 cos(x(1)) - cos(pi/(2*m))*cos(x(3))
                 sin(x(1)) - sin(pi/(2*m))/sin(x(2)/2)
                 x(4) - 2*tan(pi/(2*m))^2*cos(x(1)) ];
    [sol, fval] = fsolve(fun, [1 1 1 1], opt);
    if norm(fval) > 1e-8
        error('overcurved:shapeSolve', ...
              'fsolve did not converge for O_geom = %.6f (residual %.2e).', ...
              Og, norm(fval));
    end
    phi = real(sol(1));  omega = real(sol(2));  delta = real(sol(3));

    ysign = -1;
    if Og > O_switch
        ysign = +1;
    end

    t = linspace(-omega/2, omega/2, 800);
    c = cos(phi);  s = sin(phi);
    x_s =         (s/4).*((c-1)/(2*c+1).*sin((2*c+1).*t) ...
                        - (c+1)/(2*c-1).*sin((2*c-1).*t) - 2*sin(t));
    y_s = ysign * (s/4).*((c-1)/(2*c+1).*cos((2*c+1).*t) ...
                        + (c+1)/(2*c-1).*cos((2*c-1).*t) - 2*cos(t));
    z_s = (s*tan(phi)/4).*cos(2*t*c);

    beta = pi/2 - delta;
    x_p = x_s;
    y_p =  y_s*cos(beta) + z_s*sin(beta);
    z_p = -y_s*sin(beta) + z_s*cos(beta);

    tan_fix = tan(m/(2*pi));
    x1 = x_p;
    y1 = y_p - (y_p(1) - x_p(1)/tan_fix);
    z1 = z_p - z_p(1);

    dd     = sqrt(gradient(x1,t).^2 + gradient(y1,t).^2 + gradient(z1,t).^2);
    arclen = 2*m*trapz(t, dd);
    x = x1/arclen;  y = y1/arclen;  z = z1/arclen;

    nl    = 2*m;
    theta = pi/m;
    X = cell(nl,1);  Y = cell(nl,1);  Z = cell(nl,1);
    X{1} = x;  Y{1} = y;  Z{1} = z;
    for j = 2:nl
        X{j} = X{j-1}*cos(theta) - Y{j-1}*sin(theta);
        Y{j} = X{j-1}*sin(theta) + Y{j-1}*cos(theta);
        Z{j} = -Z{j-1};
    end
    for j = 2:nl
        X{j} = X{j} - (X{j}(1) - X{j-1}(end));
        Y{j} = Y{j} - (Y{j}(1) - Y{j-1}(end));
        Z{j} = Z{j} - (Z{j}(1) - Z{j-1}(end));
    end

    gap = norm([X{nl}(end) - X{1}(1), ...
                Y{nl}(end) - Y{1}(1), ...
                Z{nl}(end) - Z{1}(1)]);

    C = cell(nl,1);
    for j = 1:nl
        C{j} = [X{j}; Y{j}; Z{j}];
    end
end

function export_png(fig, fname, res)

    if exist('exportgraphics', 'file') == 2 || ...
       exist('exportgraphics', 'builtin') == 5
        exportgraphics(fig, fname, 'Resolution', res, ...
                       'BackgroundColor', 'white');
    else
        print(fig, fname, '-dpng', sprintf('-r%d', res));
    end
end

function p = pwd_if_relative(d)
    if isempty(d)
        p = pwd;
    elseif ~isempty(regexp(d, '^([A-Za-z]:|\\\\|/)', 'once'))
        p = d;
    else
        p = fullfile(pwd, d);
    end
end

function opt = parse_args(varargin)

    opt.Ogeom      = [];
    opt.Count      = 10;
    opt.Range      = [1.02 2.996];
    opt.OutDir     = 'shapes';
    opt.Save       = true;
    opt.View       = [45 25];
    opt.Resolution = 400;
    opt.Panel      = true;
    opt.LineWidth  = 2.2;
    opt.HighlightFirstLobe = false;

    if mod(numel(varargin), 2) ~= 0
        error('overcurved:badOptions', ...
              'Name-value options must come in pairs.');
    end

    for k = 1:2:numel(varargin)
        name = lower(char(varargin{k}));
        val  = varargin{k+1};
        switch name
            case 'ogeom',      opt.Ogeom      = double(val);
            case 'count',      opt.Count      = round(double(val));
            case 'range',      opt.Range      = double(val);
            case 'outdir',     opt.OutDir     = char(val);
            case 'save',       opt.Save       = logical(val);
            case 'view',       opt.View       = double(val);
            case 'resolution', opt.Resolution = round(double(val));
            case 'panel',      opt.Panel      = logical(val);
            case 'linewidth',  opt.LineWidth  = double(val);
            case 'highlightfirstlobe'
                               opt.HighlightFirstLobe = logical(val);
            otherwise
                error('overcurved:badOptions', ...
                      'Unknown option ''%s''.', name);
        end
    end

    if numel(opt.Range) ~= 2 || opt.Range(1) >= opt.Range(2)
        error('overcurved:badOptions', ...
              'Range must be [lo hi] with lo < hi.');
    end
    if opt.Count < 1
        error('overcurved:badOptions', 'Count must be at least 1.');
    end
    if numel(opt.View) ~= 2
        error('overcurved:badOptions', 'View must be [azimuth elevation].');
    end
end
