function varargout = ParallelRigidityTest(AdjMat,d,tolNull,mode)
% ParallelRigidityTest
%
% Polymorphic version with integrated largest parallel rigid component
%
% Usage:
%
%   isRigid = ParallelRigidityTest(A,d);
%
%       Test only. It builds the sparse parallel-rigidity matrix J, augments
%       it with the translation and scale-removal rows, and calls:
%
%           [nullity_est, s_small, stats] = estimate_nullity(Jaug,kNull,tolNull,mode)
%
%       The graph is declared parallel rigid in R^d when Jaug has zero
%       numerical nullity. This matches the logic of the original TestMat.
%
%   [isRigid,Acomp] = ParallelRigidityTest(A,d);
%
%       First performs the same test.
%       If isRigid is true, Acomp is the cleaned adjacency matrix A.
%       If isRigid is false, Acomp is the adjacency matrix induced by the
%       largest maximally parallel rigid component.
%
%   [isRigid,Acomp,testInfo] = ParallelRigidityTest(A,d);
%
%       Also returns diagnostics from estimate_nullity.
%
% Optional:
%
%   tolNull : default 1e-8
%   mode    : default 'svds'
%
% Notes:
%
%   - The test uses estimate_nullity(J,kNull,tolNull,mode).
%   - The component-extraction phase still uses gt_nulls(RTRmat), because
%     Kennedy's algorithm needs a basis of the nullspace, not only its
%     dimension.
%   - The helpers for constructing J are shared by the test and by the
%     component extraction phase.
%
% Important:
%
%   If you want to reproduce the old three-output interface
%
%       [inds, comps, Acomp] = ...
%
%   this file intentionally no longer does that. The public output is now
%   only isRigid, optionally Acomp, optionally testInfo.

if nargin < 3 || isempty(tolNull)
    tolNull = 1e-8;
end

if nargin < 4 || isempty(mode)
    mode = 'svds';
end

kNull = d+6;

% Clean adjacency matrix once.
AdjMat = clean_symmetric_adj_local(AdjMat);

% Integrated test. This is the only path used when nargout == 1.
[isRigid,testInfo] = parallel_rigidity_test(AdjMat,d,kNull,tolNull,mode);

if nargout <= 1
    varargout{1} = isRigid;
    return;
end

% If the whole graph is rigid, the largest rigid component is the whole graph.
if isRigid
    Acomp = AdjMat;
else
    % Correctly request the third output of the component routine.
    [~,~,Acomp] = LargestMaxParRigComp_component(AdjMat,d);
end

varargout{1} = isRigid;
varargout{2} = Acomp;

if nargout >= 3
    varargout{3} = testInfo;
end

end

%-------------------------------------------------------------------------------

function [isRigid,testInfo] = parallel_rigidity_test(AdjMat,d,kNull,tolNull,mode)

n = size(AdjMat,1);

[Ind_j,Ind_i] = find(tril(AdjMat,-1));
m = length(Ind_i);

testInfo = struct();
testInfo.n = n;
testInfo.m = m;
testInfo.d = d;
testInfo.kNull = kNull;
testInfo.tolNull = tolNull;
testInfo.mode = mode;
testInfo.expected_nullity = 0;

if n <= 1 || m == 0
    testInfo.nullity_est = d*n;
    testInfo.extra_nullity = d*n;
    testInfo.s_small = [];
    testInfo.stats = [];
    isRigid = false;
    return;
end

% Shared helper: same random locations and same J builder used elsewhere.
t_rand = random_centered_locations_local(d,n);
J = build_parallel_J_sparse(Ind_i,Ind_j,t_rand,d,n);

% Keep the same logic as the original test.
%
% The original code used
%
%   TestMat = J'*J + kron(ones(n,n),speye(d)) + t_rand(:)*t_rand(:)';
%   lambd   = eigs(TestMat,1,'sm');
%   isRigid = lambd >= 1e-6;
%
% Instead of forming TestMat, observe that
%   kron(ones(n,n),I_d) = C'*C,
% where
%   C = kron(ones(1,n),I_d),
% and
%   t_rand(:)*t_rand(:)' = q'*q,
% where
%   q = t_rand(:)'.
%
% Hence
%   TestMat = Jaug' * Jaug
% with
%   Jaug = [J; C; q].
%
% Therefore the original test is equivalent to checking that Jaug has
% zero nullity.
C = kron(sparse(ones(1,n)), speye(d));
q = sparse(t_rand(:).');

Jaug = [J; C; q];

[nullity_est,s_small,stats] = estimate_nullity(Jaug,kNull,tolNull,mode);

testInfo.nullity_est = nullity_est;
testInfo.s_small = s_small;
testInfo.stats = stats;

% Original-logic criterion: TestMat must be nonsingular.
isRigid = (nullity_est < 1);

end

%-------------------------------------------------------------------------------

function [LarMaxRigCompInds,MaxRigidComps,AdjLarMaxRigComp] = LargestMaxParRigComp_component(AdjMat,d)
% Component extraction part only.
%
% This is the optimized Kennedy component extraction routine. It reuses the
% same J builder used by the test; RTRmat is formed as J'*J.

tolSame = 1e-8;

n0 = size(AdjMat,1);

G0 = graph(AdjMat);
cc = conncomp(G0);
ccSizes = accumarray(cc(:),1);

MaxRigidComps = {};
bestSize = 0;
LarMaxRigCompInds = [];
AdjLarMaxRigComp = sparse([]);

[~,ordCC] = sort(ccSizes,'descend');

for ccIdx = ordCC(:).'

    nodesGlobal = find(cc == ccIdx);
    nc = numel(nodesGlobal);

    if nc <= 2
        continue;
    end

    if nc <= bestSize
        continue;
    end

    A = AdjMat(nodesGlobal,nodesGlobal);

    [localLargest,localComps] = localLargestMaxParRigComp(A,d,tolSame);

    if isempty(localLargest)
        continue;
    end

    for kk = 1:numel(localComps)
        compGlobal = nodesGlobal(localComps{kk});
        MaxRigidComps{end+1} = compGlobal(:).';
    end

    localSize = numel(localLargest);

    if localSize > bestSize
        bestSize = localSize;
        LarMaxRigCompInds = nodesGlobal(localLargest);
        LarMaxRigCompInds = LarMaxRigCompInds(:).';
        AdjLarMaxRigComp = AdjMat(LarMaxRigCompInds,LarMaxRigCompInds);
    end
end

if isempty(LarMaxRigCompInds)
    fprintf('\n No maximally rigid component of size > 2 ! \n');
    LarMaxRigCompInds = [];
    MaxRigidComps = {};
    AdjLarMaxRigComp = sparse([]);
end

end

%-------------------------------------------------------------------------------

function [LarMaxRigCompInds,MaxRigidComps] = localLargestMaxParRigComp(AdjMat,d,tolSame)

n = size(AdjMat,1);

[Ind_j,Ind_i] = find(tril(AdjMat,-1));
m = length(Ind_i);

if m == 0
    LarMaxRigCompInds = [];
    MaxRigidComps = {};
    return;
end

t_rand = random_centered_locations_local(d,n);

% Shared helper: build J once, then form RTRmat.
J = build_parallel_J_sparse(Ind_i,Ind_j,t_rand,d,n);
RTRmat = J' * J;
RTRmat = 0.5 * (RTRmat + RTRmat.');

% Use fast sparse nullspace routine.
% IMPORTANT: gt_nulls must return a nullspace basis as first output.
[Nmat,~,~] = gt_nulls(RTRmat);

if isempty(Nmat) || size(Nmat,1) ~= d*n
    error(['gt_nulls(RTRmat) must return a nullspace basis with ', ...
           'size(Nmat,1) = d*n. If it only returns nullity, ', ...
           'this algorithm cannot use it directly.']);
end

MaxRigidComps = {};
Nodes_in_max_rig_comps = [];

i = 1;
totSizeComps = 0;

while (totSizeComps < n) && (i <= n)

    indVec_i = [1:i-1,i+1:n];

    Ni = Nmat(d*(i-1)+1:d*i,:);
    Nmat_i = [Nmat(1:d*(i-1),:); Nmat(d*i+1:end,:)] ...
             - kron(ones(n-1,1),Ni);

    SameCompInds = same_component_matrix_from_Nmat_i(Nmat_i,d,tolSame);

    Nodes_i = full(sum(SameCompInds,1)) > 0;

    if any(Nodes_i)

        indVec_active = indVec_i(Nodes_i);
        S = SameCompInds(Nodes_i,Nodes_i);

        Gtmp = graph(S);
        C_Comps = conncomp(Gtmp);
        S_Comps = max(C_Comps);

        for i_k = 1:S_Comps
            comp = sort([indVec_active(C_Comps == i_k),i]);

            if numel(comp) > 2
                MaxRigidComps{end+1} = comp;
                Nodes_in_max_rig_comps = union(Nodes_in_max_rig_comps,comp);
            end
        end

        totSizeComps = length(Nodes_in_max_rig_comps);
    end

    i = i + 1;
end

if isempty(MaxRigidComps)
    LarMaxRigCompInds = [];
    return;
end

CompSizes = cellfun(@numel,MaxRigidComps);
[~,i_LarMaxComp] = max(CompSizes);
LarMaxRigCompInds = MaxRigidComps{i_LarMaxComp};

end

%-------------------------------------------------------------------------------

function J = build_parallel_J_sparse(Ind_i,Ind_j,t_rand,d,n)
% build_parallel_J_sparse
%
% Build the sparse parallel-rigidity matrix
%
%   J = blkdiag(Q_ij) * kron(B,I_d),
%
% where B is the signed edge-vertex incidence matrix and
% Q_ij = I - u_ij u_ij'.

m = length(Ind_i);

if m == 0
    J = sparse(0,d*n);
    return;
end

% Incidence-like matrix.
j_Vec_Lmat = [Ind_i Ind_j]';
j_Vec_Lmat = j_Vec_Lmat(:);

i_Vec_Lmat = kron((1:m)',ones(2,1));
val_Vec_Lmat = kron(ones(m,1),[1;-1]);

B = sparse(i_Vec_Lmat,j_Vec_Lmat,val_Vec_Lmat,m,n,2*m);
Lmat = kron(B,speye(d));

% Build block diagonal Q = blkdiag(Q_ij).
I = zeros(d*d*m,1);
Jidx = zeros(d*d*m,1);
V = zeros(d*d*m,1);

ptr = 1;
Id = eye(d);

for k = 1:m

    v = t_rand(:,Ind_i(k)) - t_rand(:,Ind_j(k));
    nv = norm(v);

    if nv < eps
        v = randn(d,1);
        nv = norm(v);
    end

    u = v / nv;
    Q = Id - u*u.';

    rows = d*(k-1) + (1:d);
    cols = rows;

    [rr,cc,val] = find(Q);

    q = numel(val);
    I(ptr:ptr+q-1) = rows(rr);
    Jidx(ptr:ptr+q-1) = cols(cc);
    V(ptr:ptr+q-1) = val;
    ptr = ptr + q;
end

I = I(1:ptr-1);
Jidx = Jidx(1:ptr-1);
V = V(1:ptr-1);

Qmat = sparse(I,Jidx,V,d*m,d*m);

J = Qmat * Lmat;

end

%-------------------------------------------------------------------------------

function SameCompInds = same_component_matrix_from_Nmat_i(Nmat_i,d,tolSame)

nminus = size(Nmat_i,1) / d;

SameCompInds = true(nminus,nminus);

for kk = 1:d

    X = Nmat_i(kk:d:end,:);

    rowNorms = sqrt(sum(X.^2,2));
    good = rowNorms > 0;

    Xn = zeros(size(X));
    Xn(good,:) = X(good,:) ./ rowNorms(good);

    C = abs(Xn * Xn.');

    SameCompInds = SameCompInds & (abs(C - 1) < tolSame);
end

SameCompInds(1:nminus+1:end) = false;
SameCompInds = sparse(SameCompInds);

end

%-------------------------------------------------------------------------------

function A = clean_symmetric_adj_local(A)

n = size(A,1);
A = spones(sparse(A));
A = spones(A | A.');
A = A - spdiags(diag(A),0,n,n);

end

%-------------------------------------------------------------------------------

function t_rand = random_centered_locations_local(d,n)

t_rand = randn(d,n);
t_rand = t_rand - mean(t_rand,2)*ones(1,n);

end




%*************************************************************************  
%  This function computes the largest maximally parallel rigid (in R^d)  
%  component of a graph with adjacency matrix AdjMat using the algorithm
%  in [1]. 
% 
%  [1] R. Kennedy, K. Daniilidis, O. Naroditsky, and C. J. Taylor, 
%      Identifying maximal rigid components in bearing-based localization, 
%      in IEEE/RSJ International Conference on Intelligent Robots and 
%      Systems, Vilamoura, Algarve, Portugal, October 2012
% 
%  Author: Onur Ozyesil
% *************************************************************************  
%  Input: 
%  AdjMat : n-by-n adjacency matrix of the graph
%  d      : dimension of the problem
%  
%  Output:
%  LarMaxRigCompInds : Indices of nodes in the largest maximally 
%                      parallel rigid component of the graph
%%*************************************************************************