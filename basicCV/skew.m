function S = skew(v)
%SKEW 3x3 skew-symmetric (cross-product) matrix: skew(u)*v = cross(u,v).
%
%   Minimal stand-in for the function of the same name in the Computer
%   Vision Toolkit (A. Fusiello), bundled so that projSFM is self-contained.
S = [   0   -v(3)  v(2);
      v(3)    0   -v(1);
     -v(2)  v(1)    0  ];
end
