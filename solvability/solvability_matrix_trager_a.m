function [S] = solvability_matrix_trager_a(G,Cams,fix_mode)
%SOLVABILITY_MATRIX_TRAGER_A Build the Trager reduced solvability matrix.
%
%   S = SOLVABILITY_MATRIX_TRAGER_A(G)
%   builds the sparse reduced finite-solvability matrix from the viewing
%   graph G using the construction from:
%
%       Federica Arrigoni, Tomas Pajdla, Andrea Fusiello,
%       "Viewing Graph Solvability in Practice", ICCV 2023.
%
%   By default, the trivial projective and scale degrees of freedom are
%   fixed by appending gauge-fixing equations.
%
%   S = SOLVABILITY_MATRIX_TRAGER_A(G, fix_mode)
%   additionally specifies the gauge-fixing mode.
%
%   INPUT:
%       G        - MATLAB graph object representing the viewing graph.
% 
%       Cams     - 3-by-4-by-n array containing camera matrices.
% 
%       fix_mode - Gauge-fixing mode:
%               'add'     append gauge/scale-fixing equations, default;
%               'remove'  remove gauge/scale columns;
%               'none'    return the raw reduced matrix.
%
%   OUTPUT:
%       S
%           Sparse reduced solvability matrix. Before gauge fixing, S has
%           size
%               11 * sum(degree(G) - 1)  by  16 * numedges(G).
%           Each edge transformation contributes 16 columns.
%
%   GAUGE CONVENTIONS:
%       In 'add' mode, the matrix is augmented with:
%           - 16 equations fixing the first edge transformation;
%           - one scale equation for each remaining edge transformation.
%
%       In 'remove' mode, the following columns are removed:
%           [1:16, 17:16:(16*numedges(G))].
%
%       This removes all 16 parameters of the first edge transformation and
%       the first parameter of every remaining edge transformation.
%
%   IMPLEMENTATION NOTES:
%       - The reduced local block has 11 rows and 16 columns.
%       - Each local block has 28 nonzero entries.
%       - Sparse triplets are preallocated and filled directly.
%
%   AUTHORS: A. Fusiello (2026), based on original implementation by 
%            F. Arrigoni (2023)


disp('Building the reduced solvability matrix...')

% Graph properties
n = numnodes(G);
m = numedges(G);
deg = degree(G);

% Recover camera centers
% C = rand(4,numnodes(G))*100;
C = zeros(4,n);
for cam = 1:n
    invV = inv3x3_fast(Cams(:,1:3,cam));
    C(:,cam) = [-invV * Cams(:,4,cam); 1];
end

% -------------------------------------------------------------------------
% Main sparse assembly
% -------------------------------------------------------------------------

nBlocks = sum(deg - 1);
nEq = 11*nBlocks;

% Each finite_equations_independent(c) has 28 nonzeros.
% Each block contributes +Ai and -Ai, hence 56 nonzeros.
nnz_per_Ai = 28;
nnzS = 2 * nnz_per_Ai * nBlocks;

I = zeros(nnzS,1);
J = zeros(nnzS,1);
V = zeros(nnzS,1);

p = 1; % current triplet index
k = 1; % current 11-row block

for vtx = 1:n

    [li,lj,lv] = finite_equations_independent_triplet(C(:,vtx));

    N = neighbors(G,vtx);
    assert(length(N) == deg(vtx))

    tau = findedge(G,vtx,N(1));

    for j = 2:deg(vtx)

        nu = findedge(G,vtx,N(j));

        row_offset_tau = 11*(k-1);
        col_offset_tau = 16*(tau-1);
        col_offset_nu = 16*(nu-1);

        q = p:p+nnz_per_Ai-1;

        I(q) = row_offset_tau + li;
        J(q) = col_offset_tau + lj;
        V(q) = lv;

        q = p+nnz_per_Ai:p+2*nnz_per_Ai-1;

        I(q) = row_offset_tau + li;
        J(q) = col_offset_nu + lj;
        V(q) = -lv;

        p = p + 2*nnz_per_Ai;
        k = k + 1;
    end
end

S = sparse(I,J,V,nEq,16*m);

% -------------------------------------------------------------------------
% Gauge fixing
% -------------------------------------------------------------------------

if nargin < 3 || isempty(fix_mode)
    fix_mode = 'add';
end

nEdgeParams = 16;
nCols = nEdgeParams * m;

switch lower(fix_mode)

    case 'none'
        % No gauge fixing.
        % S is already the raw reduced solvability matrix.

    case 'remove'
        % Remove:
        %   - all 16 columns of the first edge transformation;
        %   - the first column of every remaining edge transformation.
         colsToRemove  = [1:16, 17:16:nCols];
         S(:, colsToRemove ) = [];

    case 'add'
        % Add:
        %   - 16 equations fixing the first edge transformation;
        %   - one scale equation for each remaining edge transformation.

        B = sparse(1:16, 1:16, 1, 16, nCols);

        if m > 1
            rowsD = repelem((1:m-1)', nEdgeParams);
            colsD = nEdgeParams + (1:nEdgeParams*(m-1))';
            valsD = ones(nEdgeParams*(m-1), 1);

            D = sparse(rowsD, colsD, valsD, m-1, nCols);
        else
            D = sparse(0, nCols);
        end

        S = [S; B; D];

    otherwise
        error('Unknown fix_mode "%s". Use "add", "remove", or "none".', fix_mode);

end

end



function [i,j,v] = finite_equations_independent_triplet(c)

C1 = c(1);
C2 = c(2);
C3 = c(3);
C4 = c(4);

i = [11
    7
    6
    8
    4
    10
    3
    9
    5
    9
    2
    1
    10
    9
    11
    5
    8
    10
    2
    4
    11
    1
    7
    3
    6
    11
    9
    10];

j = [ 1
    2
    3
    4
    5
    6
    7
    7
    8
    8
    9
    10
    10
    11
    11
    12
    12
    12
    13
    13
    13
    14
    14
    15
    15
    15
    16
    16];

v = [-C4
    -C4
    -C4
    -C3
    -C4
    -C3
    -C4
    -C3
    -C3
    C4
    -C4
    -C4
    C2
    C2
    C4
    C2
    C1
    -C4
    C3
    C2
    C1
    C3
    C1
    C2
    C1
    -C3
    -C2
    C3];

end