function P = random_camera()
%RANDOM_CAMERA Generate a random full-rank camera matrix.
%
%   P = RANDOM_CAMERA() generates a random 3-by-4 camera matrix P whose
%   left 3-by-3 block is nonsingular.
%
%   Output:
%       P   - Random 3-by-4 camera matrix. The matrix is regenerated until
%             det(P(:,1:3)) has absolute value at least 1e-3.
%
%   The entries of P are sampled independently from the uniform distribution
%   on the interval (0,1).

P = rand(3,4);
while abs(det(P(1:3,1:3))) < 1e-3 ||  norm(P(:,4)) < 1e-9
    P = rand(3,4);
end
P = P/norm(P,'fro');




