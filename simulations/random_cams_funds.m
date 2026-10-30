function [Cams,Funds] = random_cams_funds(E,n,m)
%% Generate random cameras and corresponding fundamental matrices
%
%   [Cams, Funds] = random_cams_funds(E, n, m)
%
%   Generates n random full-rank 3-by-4 camera matrices and, if requested,
%   computes m noiseless fundamental matrices associated with the camera
%   pairs listed in E.
%
%   Inputs:
%       E   - m-by-2 array of camera-pair indices. Each row E(k,:) = [i j]
%             specifies that the k-th fundamental matrix is computed from
%             cameras Cams(:,:,i) and Cams(:,:,j).
%       n   - Number of cameras to generate.
%       m   - Number of fundamental matrices to generate. Typically this
%             should equal size(E,1).
%
%   Outputs:
%       Cams  - 3-by-4-by-n array containing the generated camera matrices.
%       Funds - 3-by-3-by-m array containing the fundamental matrices,
%               computed only if requested.
%
%   Cams = random_cams_funds(E,n,m); computes only the cameras, while
%   [Cams,Funds] = random_cams_funds(E,n,m); computes both

Cams = zeros(3,4,n);

for i = 1:n
    Cams(:,:,i) = random_camera();
end

if nargout > 1
    Funds = zeros(3,3,m);

    for k = 1:m
        i = E(k,1);
        j = E(k,2);
        Funds(:,:,k) = fundamental_from_P(Cams(:,:,i), Cams(:,:,j));
    end
end

end