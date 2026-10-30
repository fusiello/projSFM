function C = projective_average(M, weights, method)
%PROJECTIVE_AVERAGE Average projective matrices stored in a cell array
%
% INPUTS:
%   M       - cell array of matrices
%   weights - optional weights (vector length n)
%             default: uniform weights
%   method  - averaging method
%             'spherical' (default)
%
% OUTPUT:
%   C       - averaged matrix
%
% C = projective_average(M);                       % uniform weights
% C = projective_average(M, 'spherical');          % method only
% C = projective_average(M, w);                    % weighted
% C = projective_average(M, w, 'spherical');       % weighted + method



% -----------------------------
% Parse inputs
% -----------------------------
if nargin < 2 || isempty(weights)
    weights = [];
end

if nargin < 3 || isempty(method)
    method = 'spherical';
end

% Allow calling as projective_average(M, method)
if ischar(weights) || isstring(weights)
    method = weights;
    weights = [];
end

% -----------------------------
% Checks
% -----------------------------
if isempty(M)
    error('Input cell array M is empty.');
end

if ~iscell(M)
    error('Input M must be a cell array.');
end

n = numel(M);
sz = size(M{1});

for i = 1:n
    if ~isequal(size(M{i}), sz)
        error('All matrices must have the same size.');
    end
end

% -----------------------------
% Default weights
% -----------------------------
if isempty(weights)
    weights = ones(n,1);
else
    weights = weights(:);
    if numel(weights) ~= n
        error('weights must have length equal to numel(M).');
    end
end

% normalize weights
s = sum(weights);
if s <= 0
    weights = ones(n,1) / n;
else
    weights = weights / s;
end

% -----------------------------
% Method switch
% -----------------------------
switch lower(method)
    case {'spherical','sphere'}
        C = spherical_projective_average(M, weights);

    otherwise
        error('Unknown method: %s', method);
end
end