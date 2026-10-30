
function Ai = inv3x3_fast(A)
%INV3X3_FAST Explicit inverse of a 3x3 matrix.
%
% Uses cross products of columns.
% Assumes A is nonsingular.

a11 = A(1,1); a21 = A(2,1); a31 = A(3,1);
a12 = A(1,2); a22 = A(2,2); a32 = A(3,2);
a13 = A(1,3); a23 = A(2,3); a33 = A(3,3);

% b1 = cross(a2, a3)
b11 = a22*a33 - a32*a23;
b12 = a32*a13 - a12*a33;
b13 = a12*a23 - a22*a13;

% b2 = cross(a3, a1)
b21 = a23*a31 - a33*a21;
b22 = a33*a11 - a13*a31;
b23 = a13*a21 - a23*a11;

% b3 = cross(a1, a2)
b31 = a21*a32 - a31*a22;
b32 = a31*a12 - a11*a32;
b33 = a11*a22 - a21*a12;

detA = a11*b11 + a21*b12 + a31*b13;

s = 1 / detA;

Ai = [b11 b12 b13;
    b21 b22 b23;
    b31 b32 b33] * s;
end
