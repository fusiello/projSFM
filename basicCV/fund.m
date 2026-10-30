function F = fund(P1, P2)
%FUND Fundamental matrix from two 3x4 camera matrices (bilinear minors),
%   such that x2' * F * x1 = 0  (Hartley & Zisserman, eq. 17.3).
%
%   Minimal stand-in for the function of the same name in the Computer
%   Vision Toolkit (A. Fusiello), bundled so that projSFM is self-contained.
F = zeros(3,3);
for i = 1:3
    Xi = P1([1:i-1, i+1:3], :);
    for j = 1:3
        Yj = P2([1:j-1, j+1:3], :);
        F(j,i) = (-1)^(i+j) * det([Xi; Yj]);
    end
end
end
