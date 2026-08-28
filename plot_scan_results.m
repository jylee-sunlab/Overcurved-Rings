function fig = plot_scan_results(res, varargin)

opt = parse_args(varargin{:});

T = res.table;
Op = T.Op;
idxDep = find(T.Oeq > 1 + 1e-4, 1, 'first');

switch lower(opt.Panels)
    case 'oeq'
        layout = [1 1];
        panels = {'oeq'};
    case 'energy'
        layout = [1 1];
        panels = {'energy'};
    case 'both'
        layout = [1 2];
        panels = {'oeq', 'energy'};
    otherwise
        error('overcurved:badPanels', ...
            'Panels must be ''oeq'', ''energy'' or ''both''.');
end

fig = figure( ...
    'Color', 'w', ...
    'Name', 'Overcurved ring energy scan', ...
    'Position', [200, 200, 1000, 400]);

for k = 1:numel(panels)

    ax = subplot(layout(1), layout(2), k);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    set(ax, ...
        'TickLabelInterpreter', 'latex', ...
        'FontSize', 11);

    switch panels{k}

        case 'oeq'
            yline(ax, 1, ...
                'Color', [0.6 0.6 0.6], ...
                'LineStyle', '--');

            plot(ax, Op, T.Oeq, 'o-', ...
                'LineWidth', 1.0, ...
                'MarkerSize', 3.0, ...
                'Color', [0.00 0.00 1.00]);

            if ~isempty(idxDep)
                xline(ax, Op(idxDep), 'r--', ...
                    sprintf('$O_p = %.4f$', Op(idxDep)), ...
                    'Interpreter', 'latex', ...
                    'LabelVerticalAlignment', 'bottom', ...
                    'LabelHorizontalAlignment', 'left');
            end

            xlabel(ax, '$O_p$', 'Interpreter', 'latex');
            ylabel(ax, '$O_{\mathrm{geom}}$', 'Interpreter', 'latex');

        case 'energy'
            plot(ax, Op, T.Utotal, 'k-', 'LineWidth', 1.5);
            plot(ax, Op, T.Ub, '-', ...
                'LineWidth', 1.5, ...
                'Color', [0.15 0.35 0.75]);
            plot(ax, Op, T.Ut, '-', ...
                'LineWidth', 1.5, ...
                'Color', [0.85 0.35 0.15]);
            plot(ax, Op, T.Ua, '-', ...
                'LineWidth', 1.5, ...
                'Color', [0.20 0.55 0.30]);

            if opt.LogEnergy
                set(ax, 'YScale', 'log');
            end

            if ~isempty(idxDep)
                xline(ax, Op(idxDep), 'r--');
            end

            xlabel(ax, '$O_p$', 'Interpreter', 'latex');
            ylabel(ax, 'Energy $[\mathrm{J}]$', 'Interpreter', 'latex');
            legend(ax, ...
                {'$U$', '$U_b$', '$U_t$', '$U_a$'}, 'Interpreter', 'latex', ...
                'Location', 'best');
    end

    xlim(ax, [min(Op), max(Op)]);
    hold(ax, 'off');
end

end

function opt = parse_args(varargin)

opt.Panels = 'both';
opt.LogEnergy = true;

if mod(numel(varargin), 2) ~= 0
    error('overcurved:badOptions', ...
        'Name-value options must come in pairs.');
end

for k = 1:2:numel(varargin)
    name = lower(char(varargin{k}));
    val = varargin{k + 1};

    switch name
        case 'panels'
            opt.Panels = char(val);
        case 'logenergy'
            opt.LogEnergy = logical(val);
        otherwise
            error('overcurved:badOptions', ...
                'Unknown option ''%s''.', name);
    end
end
end
