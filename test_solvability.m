%% 
% test_solvability.m
%
% Test solvability_all on three 90-node graphs (data/A139.mat,
% data/A140.mat, data/A170.mat) for all selected methods and
% nullity-estimation modes.
%
% Expected:
%   G90-139 : summary.nullity_est > 0   (not solvable)
%   G90-140 : summary.nullity_est == 0  (solvable)
%   G90-170 : summary.nullity_est == 0  (solvable)
%
% Run setup_projSFM first (or have the toolbox on the MATLAB path).

clear; clc;

% -------------------------------------------------------------------------
% Settings
% -------------------------------------------------------------------------

rng(0);                         % reproducibility of random constructions

dataDir    = fullfile(fileparts(mfilename('fullpath')), 'data');
graphFiles = fullfile(dataDir, {'A139.mat', 'A140.mat', 'A170.mat'});
graphNames = {'G90-139',     'G90-140',     'G90-170'};

expectedNullityPositive = containers.Map();
expectedNullityPositive('G90-139')  = true;
expectedNullityPositive('G90-140')  = false;
expectedNullityPositive('G90-170') = false;

% Methods implemented in solvability_all
methods = {'pami', 'skew', 'fund'};

% Canonical modes implemented in estimate_nullity
modes = {'svd', 'svds', 'eigs_sm', 'eigs', 'auto','gt_nulls'};

% Use kSmall > 1 so that a positive nullity is less likely to be truncated.
% If the expected nullity is known to be at most 1, kSmall = 1 is enough.
kSmall = 6;

tol = 5e-8;

% -------------------------------------------------------------------------
% Run tests
% -------------------------------------------------------------------------

rows = {};
rowId = 0;

fprintf('\n============================================================\n');
fprintf('Testing solvability_all \n');
fprintf('============================================================\n\n');

for gidx = 1:numel(graphFiles)

    graphFile = graphFiles{gidx};
    graphName = graphNames{gidx};

    A = load_adjacency_from_mat(graphFile);

    fprintf('\n------------------------------------------------------------\n');
    [~, gfName, gfExt] = fileparts(graphFile);
    fprintf('Graph: %s  (%s)\n', graphName, [gfName gfExt]);
    fprintf('n = %d, m = %d\n', size(A,1), nnz(triu(A ~= 0, 1)));
    fprintf('Expected: nullity_est %s\n', ...
        ternary(expectedNullityPositive(graphName), '> 0', '== 0'));
    fprintf('------------------------------------------------------------\n');

    for midx = 1:numel(methods)

        method = methods{midx};

        for vidx = 1:numel(modes)

            mode = modes{vidx};
            rowId = rowId + 1;
            fprintf('\n[%s] method=%s, mode=%s\n', graphName, method, mode);

           % try
                summary = solvability_all(A, kSmall, tol, method, mode);
                
                nullity_est = summary.nullity_est;

                shouldBePositive = expectedNullityPositive(graphName);

                if shouldBePositive
                    pass = isfinite(nullity_est) && nullity_est > 0;
                else
                    pass = isfinite(nullity_est) && nullity_est == 0;
                end

                status = ternary(pass, 'PASS', 'FAIL');
                errMsg = "";

                fprintf('Result: nullity_est = %g --> %s\n', nullity_est, status);

            % catch ME
            % 
            %     nullity_est = NaN;
            %     pass = false;
            %     status = 'ERROR';
            %     errMsg = string(ME.message);
            % 
            %     fprintf('ERROR: %s\n', ME.message);
            % 
            % end

            rows(rowId,:) = { ...
                graphName, ...
                method, ...
                mode, ...
                nullity_est, ...
                expected_string(graphName, expectedNullityPositive), ...
                status, ...
                errMsg ...
                }; %#ok<SAGROW>

        end
    end
end

% -------------------------------------------------------------------------
% Report table
% -------------------------------------------------------------------------

Results = cell2table(rows, ...
    'VariableNames', {'graph', 'method', 'mode', 'nullity_est', ...
                      'expected', 'status', 'error_message'});

disp(' ');
disp('============================================================');
disp('Summary table');
disp('============================================================');
disp(Results);

nFail = sum(Results.status == "FAIL");
nErr  = sum(Results.status == "ERROR");

fprintf('\n============================================================\n');
fprintf('Final result\n');
fprintf('============================================================\n');
fprintf('PASS rows : %d\n', sum(Results.status == "PASS"));
fprintf('FAIL rows : %d\n', nFail);
fprintf('ERROR rows: %d\n', nErr);

if nFail == 0 && nErr == 0
    fprintf('\nALL TESTS PASSED.\n');
else
    fprintf('\nSOME TESTS FAILED OR ERRORED.\n');
end

% Optionally save results
writetable(Results, fullfile(tempdir, 'solvability_all_test_results.csv'));

% -------------------------------------------------------------------------
% Local helper functions
% -------------------------------------------------------------------------

function A = load_adjacency_from_mat(filename)
%LOAD_ADJACENCY_FROM_MAT Load an adjacency matrix from a .mat file.
%
% The .mat file may contain:
%   - an adjacency matrix A,
%   - a variable with the same name as the file, e.g. G89,
%   - a MATLAB graph object,
%   - or any single square numeric/logical matrix.

    data = load(filename);
    names = fieldnames(data);

    [~, base, ~] = fileparts(filename);

    % Prefer a variable whose name matches the file.
    if isfield(data, base)
        X = data.(base);
        A = object_to_adjacency(X);
        return;
    end

    % Prefer a variable called A.
    if isfield(data, 'A')
        X = data.A;
        A = object_to_adjacency(X);
        return;
    end

    % Otherwise search for the first usable graph or square matrix.
    for i = 1:numel(names)
        X = data.(names{i});
        try
            A = object_to_adjacency(X);
            return;
        catch
            % keep searching
        end
    end

    error('Could not find an adjacency matrix or graph object in "%s".', filename);
end

function A = object_to_adjacency(X)
%OBJECT_TO_ADJACENCY Convert a graph object or matrix to sparse adjacency.

    if isa(X, 'graph')
        A = adjacency(X);
        A = sparse(A ~= 0);
        return;
    end

    if isnumeric(X) || islogical(X)
        if ismatrix(X) && size(X,1) == size(X,2)
            A = sparse(X ~= 0);
            A = A - diag(diag(A));     % remove self-loops, if any
            A = spones(A | A.');       % force undirected pattern
            return;
        end
    end

    error('Unsupported object type for adjacency conversion.');
end

function s = expected_string(graphName, expectedNullityPositive)
    if expectedNullityPositive(graphName)
        s = "> 0";
    else
        s = "== 0";
    end
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end