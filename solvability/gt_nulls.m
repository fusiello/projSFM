function [N,flag,bound] = gt_nulls(A,k,tol)
%GT_NULLS Null space of a possibly sparse matrix.
%
%   [N,flag,bound] = GT_NULLS(A,k,tol) attempts to compute the null space
%   of A.
%
%   k
%       Initial estimate of the null-space dimension. Default: 1.
%
%   tol
%       Tolerance used to accept a vector v as null if norm(A*v) < tol.
%       Default: max(m,n)*eps(normest(A)).
%
%   N
%       Matrix whose columns are orthonormal null vectors.
%
%   flag
%       0 if the null space was computed successfully.
%       1 if the null space might be larger than the number of columns in N.
%
%   bound
%       If flag == 0, this is the dimension of null(A).
%       If flag == 1, this is an upper bound on dim(null(A)).
%
%   The algorithm uses normalized inverse iteration on U, condition
%   estimation on L1, and normalized inverse iteration on L1*U.
%
%   Reference:
%       Gotsman and Toledo,
%       "On the computation of null spaces of sparse rectangular matrices"
%       SIAM Journal on Matrix Analysis and Application 30(2):445-463
%       DOI: 10.1137/050638369
%
%   Refactored by Andrea Fusiello, 2023
%   Cleaned up version, 2026

if nargin < 1
    error('Not enough input arguments.');
end

if ~issparse(A)
    A = sparse(A);
end

[m,n] = size(A);

if m < n
    error('gt_nulls: A must be slim.');
end

if nargin < 2 || isempty(k)
    k = 1;
end

if nargin < 3 || isempty(tol)
    tol = max(m,n) * eps(normest(A));
end

warning off MATLAB:nearlySingularMatrix

t = tic;
[L,U,~,colP] = lu(A,'vector');
tfact = toc(t);

d = abs(diag(U));
d = sort(d, 'ascend');
nullity_est = n - nnz(d > tol);

fprintf(1,'Factorization time in nulls is %.2f seconds\n', tfact);
fprintf(1,'nulls: first nullity estimate = %d\n', nullity_est);

% Inverse column permutation
i_colP = zeros(1,length(colP));
i_colP(colP) = 1:length(colP);

% Replace zero diagonal entries of U by a tiny value
d = diag(U);
U = U - diag(d);
d(d == 0) = 1e-100;
U = U + diag(d);

if m > n
    L1 = L(1:n,:);
else
    L1 = L;
end

%% First stage: null space of U

kin = k;
kout = 0;
Nu = zeros(n,0);

while kout < n
    Nu = normalizedInverseTriangularIteration(U, kin, n, i_colP);
    [Nu,kout] = filterPrefixColumns(A, Nu, tol);

    if kout < kin
        break;
    end

    kin = 2 * kin;
end

dim_null_U = kout;

%% Second stage: test conditioning of L1

Nl = normalizedInverseTriangularIteration(L1, 1, n, i_colP);
[~,kout] = filterPrefixColumns(L1, Nl, tol);

if kout == 0
    % Case 1: null(A) = null(U)
    N = Nu;
    flag = 0;
    bound = dim_null_U;

    fprintf(1,'nulls: case 1, null(A)==null(U), dim(null(A))==%d\n', bound);

    warning on MATLAB:nearlySingularMatrix
    return
end

%% Third stage: null space of L1*U

fprintf(1,'nulls: L1 is ill conditioned\n');

kin = dim_null_U + 1;
kout = 0;
Nlu = zeros(n,0);

while kout < n
    Nlu = normalizedInverseLUIteration(L1, U, kin, n, i_colP);
    [Nlu,kout] = filterPrefixColumns(L1, U * Nlu, tol);

    if kout < kin
        break;
    end

    kin = 2 * kin;
end

bound = kout;

%% Final decision

if bound == dim_null_U
    % Case 2: L1 is ill conditioned, but null(A) = null(U)
    N = Nu;
    flag = 0;

    fprintf(1,'nulls: case 2, null(A)==null(U), dim(null(A))==%d\n', bound);

elseif bound == 1 && norm(A * Nlu(:,1)) > tol
    % Case 3: null(L1*U) is spurious, so null(A) is trivial
    N = zeros(n,0);
    flag = 0;
    bound = 0;

    fprintf(1,'nulls: case 3, null(A)==null(L1*U), dim(null(A))==0\n');

else
    % Cases 4 and 5: combine candidates from U and L1*U
    [N,~] = qr([Nu Nlu], 0);
    [N,kout] = filterPrefixColumns(A, N, tol);

    if kout == bound
        flag = 0;
        fprintf(1,'nulls: case 4, found null(A), dim(null(A))==%d\n', bound);
    else
        flag = 1;
        fprintf(1,'nulls: case 5, %d <= dim(null(A)) <= %d\n', kout, bound);
    end
end

warning on MATLAB:nearlySingularMatrix

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalized inverse LU iteration                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = normalizedInverseLUIteration(L,U,k,n,i_colP)

W = rand(n,k);
[W,~] = qr(W,0);

for it = 1:3
    Tt = U' \ W;
    Q  = L' \ Tt;
    [Q,~] = qr(Q,0);

    Y = L \ Q;
    X = U \ Y;
    [W,~] = qr(X,0);
end

W = sparse(W(i_colP,:));

if norm(isnan(W) + isinf(W), 1) > 0
    disp('nulls: ****************************************')
    disp('nulls: *** infs or nans in inverse iteration ***')
    disp('nulls: ****************************************')
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalized inverse iteration for a triangular Z    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = normalizedInverseTriangularIteration(Z,k,n,i_colP)

W = rand(n,k);
[W,~] = qr(W,0);

for it = 1:3
    Q = Z' \ W;
    [Q,~] = qr(Q,0);

    X = Z \ Q;
    [W,~] = qr(X,0);
end

W = sparse(W(i_colP,:));

if norm(isnan(W) + isinf(W), 1) > 0
    disp('nulls: ****************************************')
    disp('nulls: *** infs or nans in inverse iteration ***')
    disp('nulls: ****************************************')
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep only prefix columns with residual below thresh%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Wout,kout] = filterPrefixColumns(Z,W,thresh)

[~,k] = size(W);
kout = 0;

while kout < k && norm(Z * W(:,kout+1)) < thresh
    kout = kout + 1;
end

Wout = W(:,1:kout);

end