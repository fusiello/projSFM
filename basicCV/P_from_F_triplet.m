function [P1,P2,P3,a,info] = P_from_F_triplet(F12,F13,F23,varargin)
%TRIANGLE_FROM_F Recover projective cameras from a compatible triangle.
%
%   [P1,P2,P3,a,info] = P_from_F_triplet(F12,F13,F23)
%   [P1,P2,P3,a,info] = P_from_F_triplet(F12,F13,F23,tol)
%
% Input:
%   F12, F13, F23 : 3x3 fundamental matrices
%   tol           : optional numerical tolerance (default: 1e-9)
%
% Output:
%   P1, P2, P3 : 3x4 projective camera matrices
%   a          : 4x1 vector of the rank-1 correction for P3
%   info       : struct with diagnostic information
%
% Gauge:
%   P1 = [I | 0]
%   P2 = [[e21]_x F12 | e21], where e21 is the left epipole of F12
%   P3 = [[e31]_x F13 | e31] + e31 * a'
%
% Degeneracy detection:
%   The function checks whether the linear system for a has full rank 4.
%   If not, the configuration is declared degenerate. This includes the
%   collinear-camera case.
%
% Notes:
%   - the F matrices are projected to rank 2
%   - the solution for a is obtained by least squares from:
%         C + a*u' + u*a' = 0

    if nargin >= 4
        tol = varargin{1};
    else
        tol = 1e-9;
    end

    info = struct();
    info.success = false;
    info.degenerate = false;
    info.message = '';
    info.rankA = NaN;
    info.residual = NaN;
    info.sigmaA = [];
    info.collinear_suspected = false;

    % Basic size check
    if ~isequal(size(F12),[3 3]) || ~isequal(size(F13),[3 3]) || ~isequal(size(F23),[3 3])
        error('All inputs must be 3x3 matrices.');
    end

    % Reject near-zero matrices
    if norm(F12,'fro') < tol || norm(F13,'fro') < tol || norm(F23,'fro') < tol
        error('At least one input matrix is numerically zero.');
    end

    % Enforce rank-2 and normalize
    F12 = enforce_rank2(F12);
    F13 = enforce_rank2(F13);
    F23 = enforce_rank2(F23);

    F12 = F12 / norm(F12,'fro');
    F13 = F13 / norm(F13,'fro');
    F23 = F23 / norm(F23,'fro');

    % Sanity check: each F should have rank 2
    if rank(F12,tol) ~= 2 || rank(F13,tol) ~= 2 || rank(F23,tol) ~= 2
        error('At least one input matrix is not a valid rank-2 fundamental matrix.');
    end

    % Camera 1
    P1 = [eye(3), zeros(3,1)];

    % Left epipoles: e_ji satisfies e_ji'' * F_ij = 0, i.e. F_ij'' * e_ji = 0
    e21 = null(F12');
    e31 = null(F13');

    if size(e21,2) ~= 1 || size(e31,2) ~= 1
        error('Failed to compute a unique epipole from one of the fundamental matrices.');
    end

    e21 = e21(:,1);
    e31 = e31(:,1);

    if norm(e21) < tol || norm(e31) < tol
        error('Computed epipole is numerically unstable.');
    end

    e21 = e21 / norm(e21);
    e31 = e31 / norm(e31);

    % Canonical P2 and initial P3 from F12 and F13
    P2  = [skew(e21) * F12, e21];
    P30 = [skew(e31) * F13, e31];

    % Check camera rank
    if rank(P2,tol) < 3
        error('Canonical camera P2 is rank-deficient. Degenerate input.');
    end
    if rank(P30,tol) < 3
        error('Initial camera P30 is rank-deficient. Degenerate input.');
    end

    % Build the linear system for a from:
    %   C + a*u' + u*a' = 0
    C = P30' * F23 * P2 + P2' * F23' * P30;
    C = 0.5 * (C + C');   % symmetrize numerically

    u = P2' * F23' * e31; % 4x1

    % If u is (near) zero, the equation loses all dependence on a
    if norm(u) < tol
        info.degenerate = true;
        info.collinear_suspected = true;
        info.message = ['Degenerate configuration: u is numerically zero, ' ...
                        'so F23 does not constrain P3. Collinear cameras are a likely cause.'];
        P3 = NaN(3,4);
        a  = NaN(4,1);
        return;
    end

    % Unknown a = [a1;a2;a3;a4]
    % For each upper-triangular entry (i,j):
    %   C_ij + a_i*u_j + u_i*a_j = 0
    A = zeros(10,4);
    b = zeros(10,1);

    row = 1;
    for i = 1:4
        for j = i:4
            if i == j
                A(row,i) = 2*u(i);
            else
                A(row,i) = u(j);
                A(row,j) = u(i);
            end
            b(row) = -C(i,j);
            row = row + 1;
        end
    end

    s = svd(A);
    info.sigmaA = s;
    info.rankA = rank(A,tol);

    % Full rank is required for a unique solution
    if info.rankA < 4
        info.degenerate = true;
        info.collinear_suspected = true;
        info.message = ['Degenerate triangle: the linear system for a is rank deficient. ' ...
                        'This occurs in non-generic cases, including collinear camera centers.'];
        P3 = NaN(3,4);
        a  = NaN(4,1);
        return;
    end

    % Solve in least-squares sense
    a = A \ b;
    info.residual = norm(A*a - b);

    % Final P3
    P3 = P30 + e31 * a';

    % Final rank check
    if rank(P3,tol) < 3
        info.degenerate = true;
        info.message = 'Recovered P3 is rank-deficient; configuration is numerically degenerate.';
        P3 = NaN(3,4);
        a  = NaN(4,1);
        return;
    end

    info.success = true;
    info.message = 'Recovery successful.';
end

function F = enforce_rank2(F)
%ENFORCE_RANK2 Project a 3x3 matrix onto the set of rank-2 matrices.
    [U,S,V] = svd(F);
    S(3,3) = 0;
    F = U*S*V';
end
