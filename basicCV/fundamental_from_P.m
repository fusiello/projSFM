function F = fundamental_from_P(P1, P2, method)
% Compute fundamental matrix from two cameras
%
% INPUT:
%   P1, P2  - 3x4 camera matrices
%   method  - (optional) method to use:
%       'minors'        (default) uses existing fund(P1,P2)
%       'pseudo'      compute via epipoles + pseudoinverse
%
% OUTPUT:
%   F       - 3x3 fundamental matrix

    if nargin < 3 || isempty(method)
        method = 'minors';
    end

    switch lower(method)

        case 'minors'
            % elemtwise from minors (from CV toolkit)
            F = fund(P1, P2);

        case 'pseudo'
            % F = [e2]_x * P2 * pinv(P1)

            % camera center of P1 (null space)
            [~,~,V] = svd(P1);
            C1 = V(:,end);
            C1 = C1 / C1(end);

            % epipole in second image
            e2 = P2 * C1;

            % skew-symmetric matrix of e2
            e2x = [   0    -e2(3)  e2(2);
                   e2(3)     0    -e2(1);
                  -e2(2)  e2(1)     0   ];

            F = e2x * P2 * pinv(P1);

        otherwise
            error('Unknown method: %s', method);
    end

    % normalize (projective scale)
    F = F / norm(F, 'fro');
end