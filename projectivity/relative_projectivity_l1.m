function [T, AT, info] = relative_projectivity_l1(A, B, opts)
%RELATIVE_PROJECTIVITY_L1 Estimate a relative projectivity with DPCP-IRLS.
%
%   T = relative_projectivity_l1(A, B)
%   [T, AT, info] = relative_projectivity_l1(A, B, opts)
%
% Estimates an n-by-n transform such that B{k} ~ A{k} * T. A and B are
% nonempty cell arrays of equally sized real matrices (typically 3x4 cameras
% or 4x4 projectivities). Independent nonzero scales/signs are permitted.
%
% Uses the SAME linear system L as relative_projectivity.m and targets
%
%       min ||L*t||_1,   subject to ||t||_2 = 1,   T = reshape(t,n,n).
%
% DPCP-IRLS starts from the ordinary SVD solution, then repeatedly computes
% inverse-residual weights and the smallest right singular vector of the
% weighted system. Weights are per scalar equation, not per camera/view.
% With a positive delta, this minimizes a Huber-smoothed approximation:
%   rho_delta(r) = |r|                       if |r| >= delta,
%                  r^2/(2*delta) + delta/2   otherwise.
% The sphere constraint is nonconvex; global optimality is not guaranteed.
% This is an algebraic residual, not the normalized angular/geometric error.
%
% OPTIONS (all optional):
%   normalize_cameras  false  Normalize each input matrix in Frobenius norm.
%                             Set true to remove arbitrary input-scale weights.
%   fix_sign           true   Make the largest-magnitude entry of T positive.
%   max_iter           200    Maximum number of weighted SVD updates.
%   tol                1e-8   Tolerance on unit-vector change, modulo sign.
%   delta              1e-6   Residual floor, relative to ||L||_2. Smaller
%                             values approximate L1 more closely.
%
% OUTPUTS:
%   T     n-by-n transform with norm(T,'fro') == 1.
%   AT    A{k}*T, using normalized A{k} when normalize_cameras is true,
%         as in the original relative_projectivity.m.
%   info  Original diagnostics rankL, singular_values, residual (L2),
%         num_views; plus residual_l1, initial_residual_l1, iterations,
%         converged, direction_change, objective_history (actual L1),
%         smoothed_objective_history, delta (in original L residual units),
%         delta_relative, row_weights, and underdetermined.
%         Histories include initialization. Singular values refer to L.
%         Row weights are proportional to 1/max(delta,abs(L*t)), capped
%         at one. Only the smoothed objective is guaranteed to decrease
%         (up to rounding), not necessarily the actual L1 at every update.
%
% REFERENCE: Tsakiris and Vidal, Dual Principal Component Pursuit,
% JMLR 19(18), 2018, section 6.2, Algorithm 2 (codimension one).
% https://www.jmlr.org/papers/v19/17-436.html

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isstruct(opts) || ~isscalar(opts)
    error('relative_projectivity_l1:InvalidOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct('normalize_cameras', false, 'fix_sign', true, ...
    'max_iter', 200, 'tol', 1e-8, 'delta', 1e-6);
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(opts, names{j})
        opts.(names{j}) = defaults.(names{j});
    end
end
validateattributes(opts.max_iter, {'numeric'}, ...
    {'scalar','real','finite','integer','positive'}, mfilename, 'opts.max_iter');
validateattributes(opts.tol, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'opts.tol');
validateattributes(opts.delta, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'opts.delta');
validateattributes(opts.normalize_cameras, {'logical','numeric'}, ...
    {'scalar','real','finite','binary'}, mfilename, 'opts.normalize_cameras');
validateattributes(opts.fix_sign, {'logical','numeric'}, ...
    {'scalar','real','finite','binary'}, mfilename, 'opts.fix_sign');

if ~iscell(A) || ~iscell(B) || isempty(A) || numel(A) ~= numel(B)
    error('relative_projectivity_l1:InvalidInputs', ...
        'A and B must be nonempty cell arrays with the same length.');
end
nviews = numel(A);
validateattributes(A{1}, {'numeric'}, ...
    {'2d','nonempty','real','finite'}, mfilename, 'A{1}');
[m, n] = size(A{1});
dim_vec = m * n;
dim_T = n * n;
if m == 3 && n == 4 && nviews < 2
    error('relative_projectivity_l1:TooFewViews', ...
        'At least two 3x4 camera matrices are required.');
end
for k = 1:nviews
    validateattributes(A{k}, {'numeric'}, ...
        {'2d','size',[m n],'real','finite'}, mfilename, 'A{k}');
    validateattributes(B{k}, {'numeric'}, ...
        {'2d','size',[m n],'real','finite'}, mfilename, 'B{k}');
    A{k} = double(A{k});
    B{k} = double(B{k});
    na = norm(A{k}, 'fro');
    nb = norm(B{k}, 'fro');
    if na == 0 || nb == 0
        error('relative_projectivity_l1:ZeroMatrix', ...
            'Input matrices must have nonzero Frobenius norm (view %d).', k);
    end
    if opts.normalize_cameras
        A{k} = A{k} / na;
        B{k} = B{k} / nb;
    end
end

% Preserve the original projector equations and column-major vectorization.
L = zeros(dim_vec * nviews, dim_T);
Ivec = eye(dim_vec);
In = eye(n);
for k = 1:nviews
    a = B{k}(:);
    rows = (k-1)*dim_vec + (1:dim_vec);
    L(rows, :) = ((a' * a) * Ivec - a * a') * kron(In, A{k});
end
if any(~isfinite(L(:)))
    error('relative_projectivity_l1:NonfiniteSystem', ...
        'System overflow; rescale the inputs or set normalize_cameras=true.');
end

% A common scale preserves the minimizer and makes delta scale-independent.
[t, singular_values] = smallest_right_vector(L);
Lscale = singular_values(1);
if Lscale == 0
    error('relative_projectivity_l1:ZeroSystem', ...
        'The input matrices provide no nonzero projective constraints.');
end
Ln = L / Lscale;
r = Ln * t;
l1_history = zeros(opts.max_iter + 1, 1);
smooth_history = zeros(opts.max_iter + 1, 1);
l1_history(1) = sum(abs(r));
smooth_history(1) = huber_l1(r, opts.delta);
converged = false;

for iteration = 1:opts.max_iter
    % Multiplying all inverse-residual weights by delta changes no minimizer
    % and keeps the weights bounded by one, avoiding huge weighted entries.
    w = opts.delta ./ max(opts.delta, abs(r));
    Lw = bsxfun(@times, sqrt(w), Ln);
    t_next = smallest_right_vector(Lw);
    if t' * t_next < 0
        t_next = -t_next;
    end
    direction_change = norm(t_next - t);
    t = t_next;
    r = Ln * t;
    l1_history(iteration + 1) = sum(abs(r));
    smooth_history(iteration + 1) = huber_l1(r, opts.delta);
    if direction_change <= opts.tol
        converged = true;
        break;
    end
end

T = reshape(t, n, n);
T = T / norm(T, 'fro');
if opts.fix_sign
    [~, idx] = max(abs(T(:)));
    if T(idx) < 0, T = -T; end
end
AT = cell(1, nviews);
for k = 1:nviews
    AT{k} = A{k} * T;
end

if nargout > 2
    residual = L * T(:);
    info.rankL = sum(singular_values > max(size(L)) * eps(Lscale));
    info.singular_values = singular_values;
    info.residual = norm(residual);
    info.num_views = nviews;
    info.residual_l1 = norm(residual, 1);
    info.initial_residual_l1 = l1_history(1) * Lscale;
    info.iterations = iteration;
    info.converged = converged;
    info.direction_change = direction_change;
    info.objective_history = l1_history(1:iteration+1) * Lscale;
    info.smoothed_objective_history = smooth_history(1:iteration+1) * Lscale;
    info.delta = opts.delta * Lscale;
    info.delta_relative = opts.delta;
    info.row_weights = opts.delta ./ max(opts.delta, abs(Ln * T(:)));
    info.underdetermined = info.rankL < dim_T - 1;
end
end

function [v, s] = smallest_right_vector(M)
% Economy SVD omits nullspace vectors when M has fewer rows than columns.
if size(M, 1) < size(M, 2)
    [~, S, V] = svd(M);
else
    [~, S, V] = svd(M, 'econ');
end
v = V(:, end);
if nargout > 1, s = diag(S); end
end

function value = huber_l1(r, delta)
a = abs(r);
small = a < delta;
a(small) = 0.5 * (a(small).^2 / delta + delta);
value = sum(a);
end
