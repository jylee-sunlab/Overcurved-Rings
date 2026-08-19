function files = save_scan_results(res, baseName, varargin)
opt = parse_args(varargin{:});

if ~ischar(baseName) && ~isstring(baseName)
    error('overcurved:badName', 'BASENAME must be text.');
end
baseName = char(baseName);

if isempty(baseName)
    error('overcurved:badName', 'BASENAME must not be empty.');
end

folder = char(opt.Folder);
if ~isempty(folder) && ~exist(folder, 'dir')
    mkdir(folder);
end

paramTbl = parameter_table(res);
files    = {};

switch lower(opt.Format)
    case 'csv'
        files{end+1} = write_one(res.table, ...
            fullfile(folder, [baseName '_summary.csv']), opt);
        files{end+1} = write_one(paramTbl, ...
            fullfile(folder, [baseName '_parameters.csv']), opt);

        if opt.SaveDiag
            files{end+1} = write_one(struct2table_safe(res.diag), ...
                fullfile(folder, [baseName '_diagnostics.csv']), opt);
        end

        if opt.SaveFields
            names = {'s0','theta_m','t','lambda'};
            for k = 1:numel(names)
                T = field_table(res, names{k});
                files{end+1} = write_one(T, fullfile(folder, ...
                    [baseName '_field_' names{k} '.csv']), opt);
            end
        end

    case 'xlsx'
        xlsFile = fullfile(folder, [baseName '.xlsx']);
        check_overwrite(xlsFile, opt);
        if exist(xlsFile, 'file')
            delete(xlsFile);
        end

        writetable(res.table, xlsFile, 'Sheet', 'summary');
        writetable(paramTbl,  xlsFile, 'Sheet', 'parameters');

        if opt.SaveDiag
            writetable(struct2table_safe(res.diag), xlsFile, ...
                'Sheet', 'diagnostics');
        end

        if opt.SaveFields
            names = {'s0','theta_m','t','lambda'};
            for k = 1:numel(names)
                writetable(field_table(res, names{k}), xlsFile, ...
                    'Sheet', ['field_' names{k}]);
            end
        end

        files{end+1} = xlsFile;

    otherwise
        error('overcurved:badFormat', ...
              'Format must be ''csv'' or ''xlsx'' (got ''%s'').', ...
              opt.Format);
end

for k = 1:numel(files)
    fprintf('Wrote %s\n', files{k});
end
end


function opt = parse_args(varargin)
opt.Folder     = '';
opt.SaveFields = false;
opt.SaveDiag   = false;
opt.Format     = 'csv';
opt.Overwrite  = false;

if mod(numel(varargin), 2) ~= 0
    error('overcurved:badOptions', ...
          'Name-value options must come in pairs.');
end

for k = 1:2:numel(varargin)
    name = lower(char(varargin{k}));
    val  = varargin{k+1};
    switch name
        case 'folder',     opt.Folder     = val;
        case 'savefields', opt.SaveFields = logical(val);
        case 'savediag',   opt.SaveDiag   = logical(val);
        case 'format',     opt.Format     = char(val);
        case 'overwrite',  opt.Overwrite  = logical(val);
        otherwise
            error('overcurved:badOptions', ...
                  'Unknown option ''%s''.', name);
    end
end
end


function path = write_one(T, path, opt)
check_overwrite(path, opt);
writetable(T, path);
end


function check_overwrite(path, opt)
if exist(path, 'file') && ~opt.Overwrite
    error('overcurved:fileExists', ...
         ['%s already exists.  Pass ''Overwrite'', true to replace it, ' ...
          'or choose another base name.'], path);
end
end


function T = parameter_table(res)
p     = res.p;
names = fieldnames(p);
vals  = cell(numel(names), 1);

for k = 1:numel(names)
    vals{k} = num2str(p.(names{k}), '%.10g');
end

metaNames = {'codeVersion','matlabRelease','timestamp','elapsedSeconds'};
for k = 1:numel(metaNames)
    if isfield(res.meta, metaNames{k})
        v = res.meta.(metaNames{k});
        if isnumeric(v)
            v = num2str(v, '%.10g');
        end
        names{end+1} = metaNames{k};  %#ok<AGROW>
        vals{end+1}  = char(v);       %#ok<AGROW>
    end
end

T = table(names, vals, 'VariableNames', {'Parameter','Value'});
end


function T = field_table(res, name)
if ~isfield(res.fields, name)
    error('overcurved:noField', 'No field named ''%s''.', name);
end

data  = res.fields.(name);
nStep = numel(data);

first = find(~cellfun(@isempty, data), 1, 'first');
if isempty(first)
    error('overcurved:emptyField', 'Field ''%s'' is empty.', name);
end

N = numel(data{first});
M = nan(N, nStep);

for k = 1:nStep
    if ~isempty(data{k}) && numel(data{k}) == N
        M(:,k) = reshape(data{k}, [], 1);
    end
end

varNames = matlab.lang.makeValidName(compose('point_%03d', 1:nStep));

T = array2table(M, 'VariableNames', varNames);
T = addvars(T, (1:N).', 'Before', 1, 'NewVariableNames', 'gridIndex');
end


function T = struct2table_safe(s)
names = fieldnames(s);
keep  = {};
cols  = {};
n     = [];

for k = 1:numel(names)
    v = s.(names{k});
    if ~isnumeric(v) && ~islogical(v)
        continue;
    end
    v = reshape(v, [], 1);
    if isempty(n)
        n = numel(v);
    end
    if numel(v) == n
        keep{end+1} = names{k};
        cols{end+1} = double(v);
    end
end

T = table(cols{:}, 'VariableNames', keep);
end
