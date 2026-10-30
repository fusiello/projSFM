function Hn = perturb_by_angle_on_sphere(H, sigma)
% PERTURB_BY_ANGLE_ON_SPHERE Apply angular noise to a matrix/projectivity.
%
% This implements the same geometric noise model used in the Julia code:
%   1) normalize Z in Frobenius norm
%   2) sample theta = abs(N(0, sigma))
%   3) generate a random tangent direction at Z on the unit sphere
%   4) rotate Z along that tangent direction by angle theta
%
% INPUTS
%   Z      - input matrix (e.g. relative projectivity)
%   sigma  - standard deviation of angular noise, in radians
%
% OUTPUT
%   Znoisy - noisy matrix, with unit Frobenius norm
%
% Example:
%   Z = rand(4,4);
%   Z = Z / norm(Z,'fro');
%   Znoisy = addAngularNoiseProjectivity(Z, deg2rad(5));

    x = H(:);
    x = x / norm(x);

    theta = abs(sigma * randn());

    xNoisy = rotateVectorOnSphere(x, theta);

    Hn = reshape(xNoisy, size(H));
end


function xRot = rotateVectorOnSphere(x, theta)
% ROTATEVECTORONSPHERE Rotate a unit vector on the sphere by angle theta.
%
% INPUTS
%   x      - input vector
%   theta  - rotation angle in radians
%
% OUTPUT
%   xRot   - rotated unit vector
%
% Geometrically:
%   xRot = cos(theta) * x + sin(theta) * v
% where v is a random unit tangent vector orthogonal to x.

    x = x(:);
    x = x / norm(x);

    dim = length(x);

    % Generate a random tangent direction
    v = randn(dim,1);

    % Project to tangent space at x
    v = v - (x' * v) * x;

    % Retry if numerically degenerate
    nv = norm(v);
    while nv < 1e-12
        v = randn(dim,1);
        v = v - (x' * v) * x;
        nv = norm(v);
    end

    v = v / nv;

    % Move along the sphere by angular distance theta
    xRot = cos(theta) * x + sin(theta) * v;

    % Re-normalize for numerical safety
    xRot = xRot / norm(xRot);
end


