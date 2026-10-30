function [T, AT, info] = relative_projectivity(A, B, opts)
%RELATIVE_PROJECTIVITY Estimate 4x4 projective transform between two sets
%
% Goal: find T (4x4) that best aligns unit-normalized vec(A{k} * T) to unit-normalized
% vec(B{k}), with antipodal points identified.
% The two sets A and B can be cameras 3x4 or projectivities 4x4
%
%   T = relative_projectivity(A, B)
%   [T, info] = relative_projectivity(A, B, opts)
%
% INPUT:
%   A, B   cell arrays of m x n camera matrices
%            such that B{k} ~ A{k} * T
%
%   opts     struct with optional fields:
%            .normalize_cameras (default: false)
%            .fix_sign          (default: true)
%
% OUTPUT:
%   T       4x4 projective transform
%   AT      transformed A (should match B)
%   info    diagnostics (rank, singular values, residual, num_views)

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'normalize_cameras'), opts.normalize_cameras = false; end
if ~isfield(opts, 'fix_sign'), opts.fix_sign = true; end

nviews = numel(A);

if nviews ~= numel(B)
    error('A and B must have same length');
end

% infer size
[m, n] = size(A{1});
dim_vec = m * n;
dim_T   = n * n;

% --------------------
% Geometric sanity check
% --------------------
if m == 3 && n == 4 && nviews < 2
    error(['At least two 3x4 camera matrices are required to estimate ', ...
        'a 4x4 projective transformation.']);
end

% checks
for k = 1:nviews
    if ~isequal(size(A{k}), [m,n]) || ~isequal(size(B{k}), [m,n])
        error('All matrices must have same size');
    end
end
% optional normalization
if opts.normalize_cameras
    for k = 1:nviews
        A{k} = A{k} / norm(A{k}, 'fro');
        B{k} = B{k} / norm(B{k}, 'fro');
    end
end

% build system
L = zeros(dim_vec * nviews, dim_T);
Ivec = eye(dim_vec);
In   = eye(n);

for k = 1:nviews
    a = B{k}(:);

    Lk = ((a' * a) * Ivec - (a * a')) * kron(In, A{k});

    rows = (k-1)*dim_vec + (1:dim_vec);
    L(rows, :) = Lk;
end

% solve
[~, S, V] = svd(L, 'econ');
t = V(:, end);
T = reshape(t, n, n);

% normalize
T = T / norm(T, 'fro');

if opts.fix_sign
    idx = find(abs(T(:)) == max(abs(T(:))), 1);
    if T(idx) < 0
        T = -T;
    end
end


for k = 1:nviews
    AT{k} = A{k} * T ; %#ok<AGROW>
end


% diagnostics
info.rankL = rank(L);
info.singular_values = diag(S);
info.residual = norm(L * T(:));
info.num_views = nviews;
end


% function T = relative_projectivity( A, B )
% % P and PP are cell arrays of PPMs
% % T is a 4x4 projectivity s.t. B{k}  = A{k} * T
%
% %
% % for  k = 1:length(A)
% %     % same norm of a matrix made of ones
% %     A{k}  = A{k}/norm(A{k},'fro')*3.4641;
% %     B{k} = B{k}/norm(B{k},'fro')*3.4641;
% %
% % end
%
% L = [];
%
% for  k = 1:length(A)
%
%     a =  B{k}(:);
%
%     L = [L; (a'*a*eye(12)-a*a')*kron(eye(4),A{k});];
%
%
%     % if (rank(L)~=15)
%     %   warning('wrong rank!');
%     % end
%
% end
%
% [~, ~, V] = svd(L);
%
% T = reshape( V(:,end),4,4);
%
% % fprintf('relative_projectivity error: %0.5g \n',norm(maxone(P{1}*T) - maxone(PP{1}),'fro'));
%
% end
%
%
%
