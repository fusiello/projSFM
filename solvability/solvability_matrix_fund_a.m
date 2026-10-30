function S = solvability_matrix_fund_a(A,Cams,fix_mode)
%SOLVABILITY_MATRIX_FUND_A Build the fundamental-matrix solvability matrix.
%
%   S = SOLVABILITY_MATRIX_FUND_A(A)
%   builds the sparse solvability matrix associated with fundamental-matrix
%   constraints on the viewing graph whose adjacency matrix is A. Uses the
%   construction from: 
%
%       Irena Dokmanovic, Andrea Fusiello, Federica Arrigoni 
%       "Viewing-graph Solvability from the Fundamental Map", 3DV 2027.
% 
%   By default, the projective gauge is fixed by removing selected columns.
%
%   S = SOLVABILITY_MATRIX_FUND_A(A, fix_mode)
%   additionally specifies the gauge-fixing mode.
%
%   INPUT:
%       A
%           n-by-n adjacency matrix of an undirected viewing graph. Only the
%           sparsity pattern is used.
% 
%       Cams     - 3-by-4-by-n array containing camera matrices.
% 
%       fix_mode
%           Gauge-fixing mode:
%               'remove'  remove gauge columns, default;
%               'add'     append gauge-fixing equations;
%               'none'    return the raw, ungauged matrix.
%
%   OUTPUT:
%       S
%           Sparse solvability matrix. Before gauge fixing, S has size
%               8*m by 11*n,
%           where n is the number of cameras and m is the number of
%           undirected edges.
%
%   CONVENTIONS:
%       Each camera contributes 11 columns, corresponding to the first
%       11 entries of vec(P), with MATLAB column-major vectorization.
%       The twelfth camera entry is removed by the local affine chart.
%
%       The gauge columns are
%
%           [1:11, 12, 15, 18, 21],
%
%       namely all 11 parameters of the first camera plus the first row of
%       the second camera, under column-major vectorization.
%
%   IMPLEMENTATION NOTES:
%       - Per-camera inverses V^{-1} and centers c = -V^{-1}m are
%         precomputed once per camera.
%       - Sparse triplets are preallocated and filled edge by edge.
%       - The local Jacobian routine returns only the retained 8-by-11
%         blocks.
%
%   ALGORITHM: Jacobian of the fundamental map (Irena's thesis)
%
%   AUTHORS: A. Fusiello (2026), based on original implementation by 
%            Irena Dokmanovic (2026)


disp('Building the solvability matrix (fund) ...')

% Graph properties
B = adj2inc(A, false);
[n, m] = size(B);
E = inc2edgelist(B);

% -------------------------------------------------------------------------
% Precompute per-camera quantities
% -------------------------------------------------------------------------
invVs   = zeros(3,3,n);
centers = zeros(3,n);

for cam = 1:n
    V = Cams(:,1:3,cam);
    mm = Cams(:,4,cam);
    invVs(:,:,cam) = inv3x3_fast(V);
    centers(:,cam) = -invVs(:,:,cam) * mm;
end

% -------------------------------------------------------------------------
% Main sparse assembly
% -------------------------------------------------------------------------
entriesPerBlock = 8 * 11;
entriesPerEdge  = 2 * entriesPerBlock;
nnzMax = entriesPerEdge * m;

rows    = zeros(nnzMax, 1);
columns = zeros(nnzMax, 1);
values  = zeros(nnzMax, 1);

% Pattern for one 8 x 11 block, in the same order as J.'
localRows = repelem((1:8)', 11);
localCols = repmat((1:11)', 8, 1);

ptr = 1;

for k = 1:m
    i = E(k, 1);
    j = E(k, 2);

    [J1, J2] = fundamental_jacobian_optim(...
        invVs(:,:,i), invVs(:,:,j), centers(:,i), centers(:,j));

    edgeRows = 8 * (k - 1) + localRows;

    idx1 = ptr : ptr + entriesPerBlock - 1;
    idx2 = ptr + entriesPerBlock : ptr + entriesPerEdge - 1;

    rows(idx1)    = edgeRows;
    columns(idx1) = 11 * (i - 1) + localCols;
    values(idx1)  = reshape(J1.', [], 1);

    rows(idx2)    = edgeRows;
    columns(idx2) = 11 * (j - 1) + localCols;
    values(idx2)  = reshape(J2.', [], 1);

    ptr = ptr + entriesPerEdge;
end

S = sparse(rows, columns, values, 8*m, 11*n);


% -------------------------------------------------------------------------
% Gauge fixing
% -------------------------------------------------------------------------

if nargin < 3 || isempty(fix_mode)
    fix_mode = 'remove';
end

% Columns used to fix the projective gauge.
gaugeCols = [1:11 12 15 18 21]; % first camera + first row of the second camera

switch lower(fix_mode)

    case 'none'
        % Do nothing.

    case 'remove'
        % Fix the gauge by removing selected columns.
        S(:, gaugeCols) = [];

    case 'add'
        % Fix the gauge by adding linear equations x(gaugeCols)=0.
        nCols = size(S, 2);
        q = numel(gaugeCols);
        G = sparse(1:q, gaugeCols, 1, q, nCols);
        S = [S; G];

    otherwise
        error('Unknown fix_mode "%s". Use "remove", "add", or "none".', fix_mode);

end

end


% =========================================================================
% Local optimized Jacobian: returns only first 11 columns
% =========================================================================

function [JF_P1, JF_P2] = fundamental_jacobian_optim(invV1,invV2,c1,c2)
%FUNDAMENTAL_JACOBIAN_OPTIM Optimized local Jacobian for F constraints.
%
%   [JF_P1, JF_P2] = FUNDAMENTAL_JACOBIAN_OPTIM(invV1, invV2, c1, c2)
%
%   computes the two 8-by-11 Jacobian blocks of the normalized fundamental
%   matrix constraint with respect to two cameras
%
%       P1 = [V1, m1],    P2 = [V2, m2],
%
%   using precomputed inverses and camera centers:
%
%       invV1 = V1^{-1},
%       invV2 = V2^{-1},
%       c1    = -V1^{-1} m1,
%       c2    = -V2^{-1} m2.
%
%   The fundamental matrix is represented as
%
%       F = V2^{-T} [c1 - c2]_x V1^{-1}.
%
%   The output uses the affine chart obtained by normalizing by F(3,3)
%   and dropping the ninth entry of vec(F). Only the first 11 camera
%   parameters are returned for each camera, matching the 11-column
%   per-camera layout used in SOLVABILITY_MATRIX_FUND_A.
%
%   INPUT:
%       invV1, invV2
%           Precomputed inverses of V1 and V2.
%
%       c1, c2
%           Precomputed camera centers in the affine chart:
%           c = -V^{-1} m.
%
%   OUTPUT:
%       JF_P1
%           8-by-11 Jacobian block with respect to camera 1.
%
%       JF_P2
%           8-by-11 Jacobian block with respect to camera 2.
%
%   Notes:
%       - This routine avoids explicit inv calls.
%       - Small Kronecker products are replaced by specialized 3-by-3
%         block constructions.
%       - It assumes V1 and V2 are nonsingular.

% ---------------------------------------------------------------------
% Basic matrices
% ---------------------------------------------------------------------

% invV1 = inv3x3_fast(V1);
% invV2 = inv3x3_fast(V2);
% c1 = -invV1*m1;
% c2 = -invV2*m2;

% For any matrix M: M*K33 == M(:,Kcols33)
% Thus multiplication by K33 can be replaced by column indexing.
% [~,Kcols33] = max(K33,[],1);
Kcols33 = [1  4  7  2  5  8  3  6  9];

d = c1 - c2;

invV2t = invV2.';

Dx = [  0    -d(3)   d(2);
    d(3)  0     -d(1);
    -d(2)  d(1)   0    ];

TS = kron_AB_times_skew_lift(invV1, invV2t);

% ---------------------------------------------------------------------
% P1 Jacobian
% ---------------------------------------------------------------------
invV2t_Dx_invV1 = invV2t * Dx * invV1;

Kc1_invV1 = [c1(1)*invV1, c1(2)*invV1, c1(3)*invV1];
J_V1 = -TS * Kc1_invV1 -kron_At_B_3x3(invV1, invV2t_Dx_invV1);
J_m1 = -TS * invV1;
JF_P1_raw = [J_V1, J_m1];

% ---------------------------------------------------------------------
% P2 Jacobian
% ---------------------------------------------------------------------
A2 = invV1.' * Dx.' * invV2;

M = kron_A_B_3x3(A2, invV2t);
Kc2_invV2 = [c2(1)*invV2, c2(2)*invV2, c2(3)*invV2];
J_V2 = -M(:, Kcols33) + TS * Kc2_invV2;
J_m2 = TS * invV2;
JF_P2_raw = [J_V2, J_m2];

% ---------------------------------------------------------------------
% Normalization and affine chart
% ---------------------------------------------------------------------
F = invV2t_Dx_invV1;

f9 = F(9); f8 = F(1:8)';
alpha = 1 / f9;
beta  = 1 / (f9 * f9);

% Since e33.'*J selects row 9 of J, and row 9 is discarded afterward,
% only the retained rows are computed.
JF_P1 = alpha * JF_P1_raw(1:8,1:11) - beta * f8 * JF_P1_raw(9,1:11);
JF_P2 = alpha * JF_P2_raw(1:8,1:11) - beta * f8 * JF_P2_raw(9,1:11);

end



function TS = kron_AB_times_skew_lift(invV1, invV2t)
% Computes:
%   TS = kron(invV1.', invV2t) * S
% without forming kron(invV1.', invV2t).
%
% S*x = vec(skew(x)), MATLAB column-major.

E1 = [ 0  0  0;
    0  0 -1;
    0  1  0];

E2 = [ 0  0  1;
    0  0  0;
    -1  0  0];

E3 = [ 0 -1  0;
    1  0  0;
    0  0  0];

TS = [ ...
    reshape(invV2t * E1 * invV1, [], 1), ...
    reshape(invV2t * E2 * invV1, [], 1), ...
    reshape(invV2t * E3 * invV1, [], 1) ...
    ];
end


function K = kron_At_B_3x3(A, B)
% Computes kron(A.', B) for 3x3 matrices A,B.
% Equivalent to kron(A.', B), but avoids kron overhead.

K = [ ...
    A(1,1)*B, A(2,1)*B, A(3,1)*B; ...
    A(1,2)*B, A(2,2)*B, A(3,2)*B; ...
    A(1,3)*B, A(2,3)*B, A(3,3)*B ...
    ];
end

function K = kron_A_B_3x3(A, B)
% Computes kron(A,B) for 3x3 matrices A,B.

K = [ ...
    A(1,1)*B, A(1,2)*B, A(1,3)*B; ...
    A(2,1)*B, A(2,2)*B, A(2,3)*B; ...
    A(3,1)*B, A(3,2)*B, A(3,3)*B ...
    ];
end