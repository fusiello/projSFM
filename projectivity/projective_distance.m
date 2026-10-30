function err = projective_distance(Hest, Htrue, method)
%PROJECTIVE_DISTANCE Compare projective matrices up to global scale/sign
%
%   err = projective_distance(Hest, Htrue)
%   err = projective_distance(Hest, Htrue, method)
%
% INPUT
%   Hest, Htrue : either
%       - single matrices
%       - cell arrays of matrices
%
%   method      : distance type
%       'sphere'  - spherical [rad] distance between normalized vectorizations       
%       'op'      - orthogonal projector distance
%       'fro'     - normalized Frobenius distance with optimal scalar alignment
%
% OUTPUT
%   err         : scalar if Hest,Htrue are single matrices
%                 column vector if they are cell arrays
%                 eith 'sphere' it is angle in rad

    if nargin < 3 || isempty(method)
        method = 'sphere';
    end

    % Wrap single matrices into cell arrays
    if ~iscell(Hest)
        Hest = {Hest};
    end

    if ~iscell(Htrue)
        Htrue = {Htrue};
    end

    n = numel(Hest);

    if n ~= numel(Htrue)
        error('Hest and Htrue must have the same length.');
    end

    err = zeros(n,1);

    for i = 1:n
        A = Hest{i};
        B = Htrue{i};

        if ~ismatrix(A) || ~ismatrix(B)
            error('Entries of Hest and Htrue must be matrices.');
        end

        if ~isequal(size(A), size(B))
            error('Size mismatch at index %d.', i);
        end

        err(i) = projective_distance_single(A, B, method);
    end

    % Return scalar when inputs were single matrices
    if n == 1
        err = err(1);
    end
end


function d = projective_distance_single(A, B, method)
%PROJECTIVE_DISTANCE_SINGLE Compare two projective matrices up to scale/sign
%
%   A and B can be any size, e.g.
%       F : 3x3
%       P : 3x4
%       T : 4x4
%
%   method options:
%       'sphere'  - distance on the unit sphere [rad]
%       'op'      - orthogonal projector distance
%       'fro'     - Frobenius distance with optimal scalar alignment


    a = A(:);
    b = B(:);

    switch lower(method)
        case {'sphere','spherical'}
            d = spherical_vec_distance(a, b);

        case 'op'
            d = op_vec_distance(a, b);

        case 'fro'
            denom = sum(a .* a);
            if denom == 0
                error('First input matrix has zero norm.');
            end
            nb = norm(B, 'fro');
            if nb == 0
                error('Second input matrix has zero norm.');
            end

            alpha = sum(a .* b) / denom;
            d = norm(alpha * A - B, 'fro') / nb;

        otherwise
            error('Unknown method: %s', method);
    end
end



function c = op_vec_distance(a , b)

%  orthogonal projection distance

% normalize and vectorize
a = a(:)/norm(a(:));
b = b(:)/norm(b(:));

c =  norm((eye(length(a)) - a * a') * b) ;
end


function c = spherical_vec_distance(a, b)

% geodesic distance with antipodal points identifyed
% it can be interpreted as an angle in rad

% normalize and vectorize
a = a(:)/norm(a(:));
b = b(:)/norm(b(:));

c = real(acos(abs(a'*b)));

end