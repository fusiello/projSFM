function S = solvability_matrix_skew_a(A, Cams, fix_mode)
%SOLVABILITY_MATRIX_SKEW_A Build the skew-symmetry solvability matrix.
%
%   S = SOLVABILITY_MATRIX_SKEW_A(A)
%   builds the sparse finite-solvability matrix associated with the
%   all-edges skew-symmetry constraints of the viewing graph whose adjacency
%   matrix is A. Uses the construction from:
%
%    Federica Arrigoni, Kathlen Kohn, Andrea Fusiello, Tomas Pajdla,
%    "An Algebraic Geometry Approach to Viewing Graph Solvability", PAMI 2026.
% 
%   By default, the projective gauge and per-camera scale
%   ambiguities are fixed by appending linear equations.
%
%   S = SOLVABILITY_MATRIX_SKEW_A(A, fix_mode)
%   additionally specifies the gauge-fixing mode.
%
%   INPUT:
%       A
%           n-by-n adjacency matrix of an undirected viewing graph. Only the
%     
%       Cams     - 3-by-4-by-n array containing camera matrices.
% 
%       fix_mode
%           Gauge-fixing mode:
%               'add'     append gauge/scale-fixing equations, default;
%               'remove'  remove gauge/scale columns;
%               'none'    return the raw, ungauged matrix.
%
%   OUTPUT:
%       S
%           Sparse solvability matrix. Before gauge fixing, S has size
%               10*m by 12*n,
%           where n is the number of cameras and m is the number of
%           undirected edges.
%
%   CONSTRAINTS:
%       The edge constraints are obtained by linearizing the skew-symmetry
%       equations associated with
%           P_j' * F_ij * P_i.
%       Each incident edge-camera pair contributes one 10-by-12 local block.
%       The implementation stores only the known nonzero entries of these
%       blocks.
%
%   GAUGE CONVENTIONS:
%       Each camera contributes 12 columns, corresponding to vec(P) for a
%       3-by-4 camera matrix P, using MATLAB column-major vectorization.
%
%       The projective gauge columns are
%           [1:12, 13, 16, 19, 22],
%       namely all 12 parameters of the first camera plus the first row of
%       the second camera. Under column-major vectorization, the first row
%       of a 3-by-4 camera has local indices [1 4 7 10].
%
%       In 'add' mode, the matrix is augmented with:
%           - 16 equations fixing the projective gauge variables;
%           - one scale equation for each camera except the first one.
%
%       In 'remove' mode, the gauge columns and one scale column per camera
%       are removed. Since column 12 belongs both to the projective gauge
%       and to the per-camera scale columns, duplicate column indices are
%       removed with UNIQUE.
%
%   IMPLEMENTATION NOTES:
%       - Sparse triplets are preallocated and filled directly.
%       - The fixed 10-by-12 block sparsity pattern is used to avoid storing
%         known zero entries.
%       - Small Kronecker products with identity/ones are replaced by local
%         helper routines.
%
%   AUTHORS: A. Fusiello (2026), based on original implementation by 
%            F. Arrigoni and A. Fusiello (2024)

disp('Building the solvability matrix (skew faster)...')


% Graph properties
B = adj2inc(A, false);
n =  size(B,1);
E = inc2edgelist(B);

% recover fundamental matrices
Funds = funds_from_cams(Cams,E);

% -------------------------------------------------------------------------
% Main sparse assembly
% -------------------------------------------------------------------------

% Internal convention: B is interpreted as edges x nodes.
B = B';

% Find the nonzero positions in the incidence matrix B
[non_zero_rows, non_zero_cols] = find(B);

% Initialize the new matrix with larger dimensions
new_rows = size(B, 1) * 10;
new_cols = size(B, 2) * 12;
% B_grande_full = zeros(nuove_righe, nuove_colonne);
% block_size = 10 * 12;

block_pattern = [1 5 6 7 11 15 16 17 21 25 26 27 32 35 38 39 ...
    42 45 48 49 52 55 58 59 63 66 68 70 73 76 78 80 ...
    83 86 88 90 94 97 99 100 104 107 109 110 114 117 119 120];

block_size_reduced = length(block_pattern);

N = nnz(B)*block_size_reduced;
rows = zeros(N,1);
cols = zeros(N,1);
val = zeros(N,1);

K = (commutation(4,4) + eye(16));

% Replace the 10x12 blocks in the new matrix
for k = 1:length(non_zero_rows)

    riga_inizio_fine =   ((non_zero_rows(k) - 1) * 10 + 1 : non_zero_rows(k) * 10)';
    colonna_inizio_fine = ((non_zero_cols(k) - 1) * 12 + 1 :  non_zero_cols(k) * 12)';

    edge = non_zero_rows(k);
    node = non_zero_cols(k); % this is one end of the edge

    % other_node = find( B(edge,:) == -B(edge, node));  % other end of the edge

    nn = find(B(edge,:));
    if node == nn(1)
        other_node = nn(2);
    else
        other_node = nn(1);
    end

    % Compute the 10x12 matrix
    F = Funds(:,:,edge); % fundamental matrix
    if B(edge, node) >   0
        Q = Cams(:,:,other_node); % left camera
        %DP = K * kron(eye(4), (F*Q)');
        DP = K*fast_kron_eye_A(4,(F*Q)');
        blocco = DP([1,6,11,16,2,3,4,7,8,12],:);
    else
        P = Cams(:,:,other_node); % right camera
        %DQ = K * kron(eye(4), P'*F);
        DQ = K*fast_kron_eye_A(4,P'*F);
        blocco = DQ([1,6,11,16,2,3,4,7,8,12],:);
    end

    % Replace the block in the new matrix
    % B_grande_full(riga_inizio_fine,  colonna_inizio_fine) = blocco;

    %r = kron(ones(12,1), (riga_inizio_fine));
    r = fast_kron_ones_A(12,1,riga_inizio_fine);

    %c = kron((colonna_inizio_fine), ones(10,1));
    c = fast_kron_A_ones(colonna_inizio_fine, 10, 1);

    rows((k - 1) * block_size_reduced + 1 : k * block_size_reduced) = r(block_pattern);
    cols((k - 1) * block_size_reduced + 1 : k * block_size_reduced) =   c(block_pattern);
    val( (k - 1) * block_size_reduced + 1 : k * block_size_reduced) = blocco(block_pattern);

end

S = sparse(rows, cols, val, new_rows, new_cols);


% -------------------------------------------------------------------------
% Gauge / scale fixing
% -------------------------------------------------------------------------

if nargin < 3 || isempty(fix_mode)
    fix_mode = 'add';
end

% Gauge:
%   - all 12 parameters of camera 1
%   - first row of camera 2, assuming MATLAB column-major vec(P)
%
% First row local indices of a 3x4 camera under column-major vectorization:
% correspond to: [1 4 7 10]
gaugeCols = [1:12, 13, 16, 19, 22];

% One scale variable per camera: remove/fix the last entry of each camera.
scaleCols = 12:12:(12*n);

switch lower(fix_mode)

    case 'none'
        % Do nothing.

    case 'remove'
        % Remove all gauge and scale columns at once.

        colsToRemove = unique([gaugeCols, scaleCols]);
        % Use unique() because column 12 belongs both to the gauge and to
        % the per-camera scale columns.
        S(:, colsToRemove) = [];

    case 'add'
        % Add linear equations fixing:
        %   - the 16 projective gauge variables;
        %   - one scale equation per camera except camera 1, because camera 1
        %     is already fully fixed by the gauge.
        nCols = size(S, 2);

        % Projective gauge equations.
        dof = numel(gaugeCols); % 16
        G = sparse(1:dof, gaugeCols, 1, dof, nCols);

        % Scale equations: one per camera, except the first camera.
        % because it is already fixed.
        D = kron(speye(n-1), sparse(ones(1,12)));
        D = [sparse(n-1,12), D];

        S = [S; G; D];

    otherwise
        error('Unknown fix_mode "%s". Use "remove", "add", or "none".', fix_mode);
end


end


function B = fast_kron_eye_A(n,A)
% B = kron(eye(n), A);
B = zeros([size(A) n n]);
B(:,:,1:n+1:n^2) = repmat(A, [1 1 n]);
B = permute(B, [1 3 2 4]);
B = reshape(B, n*size(A));
end

function B = fast_kron_ones_A(m,n,A)
% B = kron(ones(m,n), A);
B = repmat(A, [m n]);
end

function B = fast_kron_A_ones(A, m, n)
% B = kron(A, ones(m,n));
p = size(A,1);
q=size(A,2);
i = 1:p;  i = i(ones(1,m),:);
j = 1:q;  j = j(ones(1,n),:);
B = A(i,j);
end

