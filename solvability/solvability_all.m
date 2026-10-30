function summary = solvability_all(A, kSmall, tol, method, mode)
%SOLVABILITY_ALL Build a solvability matrix and estimate its numerical nullity.
%
%   summary = SOLVABILITY_ALL(A, kSmall, tol, method, mode)
%   constructs one of several solvability matrices associated with the
%   input graph A, row-normalizes the resulting matrix, and estimates its
%   numerical nullity from its smallest singular values.
%
%   The graph is represented by its adjacency matrix A. The function first
%   converts A to an incidence matrix and an edge list, then builds a
%   method-dependent solvability matrix S. The matrix S is row-normalized
%   before calling ESTIMATE_NULLITY.
%
%   Inputs:
%       A
%           Adjacency matrix of the graph. The graph is treated as
%           undirected by the call
%               B = adj2inc(A, false);
%
%       kSmall
%           Number of smallest singular values to compute or estimate.
%           Default: 1.
%
%       tol
%           Numerical threshold for declaring a singular value to be zero.
%           The estimated nullity is
%               #{ i : s_small(i) < tol }.
%           Default: 5e-8.
%
%       mode
%           Singular-value estimation mode passed to ESTIMATE_NULLITY.
%           Typical values are:
%
%               'svd'
%                   Full dense economy SVD.
%
%               'svds'
%                   Sparse computation of the smallest singular values.
%
%               'eigs'
%                   Eigenvalue-based estimate using J'*J and shift-invert.
%
%               'eigs_sm'
%                   Eigenvalue-based estimate using smallest-magnitude
%                   eigenvalues of J'*J.
%
%               'auto'
%                   Try the default eigs strategy first, then fall back to
%                   shift-invert.
%
%               'gt_nulls'
%                   Null-space dimension via GT_NULLS (Gotsman-Toledo).
%
%           Default: depends on method ('svds' for trager, 'eigs' for
%           4dof, 'eigs_sm' otherwise).
%
%       NOTE: method comes BEFORE mode in the argument list.
%
%       method
%           Solvability-matrix construction method. Supported values are:
%
%               'trager'
%                   ICCV23 reduced Trager-Osserman-Ponce construction,
%                   using SOLVABILITY_MATRIX_TRAGER_A.
%
%               'pami'
%                   Legacy ECCV24/PAMI25 all-edges skew-symmetry
%                   construction, using SOLVABILITY_MATRIX_SKEW.
%                   The random sampling is internal to that routine.
%
%               'skew'
%                   Optimized all-edges skew-symmetry construction, using
%                   RANDOM_CAMS_FUNDS and SOLVABILITY_MATRIX_SKEW_A.
%
%               'fund'
%                   Optimized fundamental-matrix construction, using
%                   RANDOM_CAMS_FUNDS and SOLVABILITY_MATRIX_FUND_A.
%
%               '4dof'
%                   Four-degree-of-freedom construction, using
%                   SOLVABILITY_MATRIX_4DOF.
%
%           Default: 'skew'.
%
%   Output:
%       summary
%           Structure containing the result and diagnostic information.
%           Mainly: summary.nullity_est
%
%   Solvability criterion:
%       The graph is considered numerically solvable when
%
%           summary.nullity_est < 1.
%
%       Equivalently, none of the computed smallest singular values lies
%       below the threshold tol. This conclusion is only as reliable as the
%       chosen value of kSmall: if kSmall is too small, the estimated nullity
%       may be truncated.
%
%   Defaults:
%       kSmall = 1
%       tol    = 5e-8
%       method = 'skew'
%       mode   = method-dependent (see above)
%
%   See also ESTIMATE_NULLITY, ADJ2INC, INC2EDGELIST, RANDOM_CAMS_FUNDS,
%   SOLVABILITY_MATRIX_SKEW, SOLVABILITY_MATRIX_SKEW_A, SOLVABILITY_MATRIX_TRAGER_A,
%   SOLVABILITY_MATRIX_FUND_A, SOLVABILITY_MATRIX_4DOF.


% these are prefererd method-mode pairs
defaultModeByMethod = {
'trager', 'svds'    
'pami',   'eigs_sm'
'skew',   'eigs_sm'
'fund',   'eigs_sm'
 '4dof',  'eigs'
    };

if nargin < 2 || isempty(kSmall), kSmall = 1; end
if nargin < 3 || isempty(tol),    tol = 5e-8; end
if nargin < 5 || isempty(method), method = 'skew'; end
if nargin < 4 || isempty(mode)
    idx = strcmpi(method, defaultModeByMethod(:,1));
    if any(idx)
        mode = defaultModeByMethod{idx,2};
    else
        mode = 'eigs_sm';
    end
end


summary = struct();
summary.method       = method;
summary.mode         = mode;

fprintf('Testing solvability on the %s constraints with %s\n', summary.method, summary.mode)

% random sampling cameras
Cams = random_cams_funds([], size(A,1));



switch method

    case{'trager'} % Trager–Osserman–Ponce construction
        % Our ICCV 23 reduced version of Trager–Osserman–Ponce
        G=graph(A);
        S = solvability_matrix_trager_a(G,Cams);

    case {'pami'} % skew implicit constraints
        % The original ECCV24/PAMI26 published version (legacy code)
        G=graph(A);
        [S,~,~,~] = solvability_matrix_skew(G,Cams);

    case {'skew'} % skew implicit constraints
         % The ECCV24/PAMI26 construction optimized by Andrea
        S = solvability_matrix_skew_a(A,Cams);

    case{'fund'} % fundamental mapp
        % The 3DV Irena's construction optimized by Andrea
        S = solvability_matrix_fund_a(A,Cams);

    case{'4dof'} % skew implicit constraints + parametrization
        % New version of skew-constraints with tree-based 4-dof 
        % per camera parametrization  
        S = solvability_matrix_4dof(A);


    otherwise
        error('method not implemented')
end

% row normalization with guard against zero rows
nr = vecnorm(S,2,2); nr(nr==0) = 1;
S = S ./ nr;

[eq,vars] = size(S);
summary.J            = S;
summary.n            = size(A,1);
summary.m_undirected = nnz(triu(A,1));
summary.eq           = eq;
summary.vars         = vars;

fprintf('n = %d\n', summary.n);
fprintf('undirected edges =  %d\n', summary.m_undirected);
fprintf('vars = %d\n', summary.vars);


%check nullity
[nullity_est, s_small, stats] = estimate_nullity(S,kSmall,tol,mode);

summary.nullity_est  = nullity_est;
summary.s_small      = s_small;
summary.thr          = stats.thr;
summary.kUse         = stats.kUse;
summary.tol          = tol;
summary.kSmall       = kSmall;
summary.smax         = stats.smax;


fprintf('estimated nullity (#small sing < %.1e) = %d (using k=%d)\n',...
    summary.thr, summary.nullity_est, summary.kUse);
disp('----------------------------------------------------');

end

