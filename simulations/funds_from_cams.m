function Funds = funds_from_cams(Cams,E,m)
%% Compute fundamental matrices from given cameras
%
%   Funds = funds_from_cams(Cams,E,m)
%
%   Computes m noiseless fundamental matrices associated with the camera
%   pairs listed in E.
%
%   Inputs:
%       Cams - 3-by-4-by-n array containing camera matrices.
%       E    - m-by-2 array of camera-pair indices. Each row E(k,:) = [i j]
%              specifies that the k-th fundamental matrix is computed from
%              cameras Cams(:,:,i) and Cams(:,:,j).
%       m    - Number of fundamental matrices to compute. Typically this
%              should equal size(E,1).
%
%   Output:
%       Funds - 3-by-3-by-m array containing the fundamental matrices.
%
%   See also fundamental_from_P.

if nargin < 3
    m = size(E,1);
end

Funds = zeros(3,3,m);

for k = 1:m
    i = E(k,1);
    j = E(k,2);

    Funds(:,:,k) = fundamental_from_P(Cams(:,:,i), Cams(:,:,j));
end

end