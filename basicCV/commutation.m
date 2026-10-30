function K = commutation(m, n)
%COMMUTATION Commutation matrix K (mn x mn) such that K*vec(A) = vec(A.')
%   for any m x n matrix A.
%
%   Minimal stand-in for the function of the same name in the Computer
%   Vision Toolkit (A. Fusiello), bundled so that projSFM is self-contained.
if nargin < 2, n = m; end
I = reshape(1:m*n, m, n);
I = I.';
K = eye(m*n);
K = K(I(:), :);
end
