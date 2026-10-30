
% test with cameras
n = 3; m = 3;
A = rand(n,m);
B = rand(n,m);

% invarianza a scala con segno
projective_distance(A(:), B(:),'sphere')
projective_distance(A(:), -10*B(:),'sphere')

projective_distance(A(:), B(:),'op')
projective_distance(A(:), -10*B(:),'op')

% invarianza a scala con segno
projective_average({A,B}, 'sphere')
projective_average({A,-10*B}, 'sphere')

% % media di A e -A deve fare A
projective_distance(A, projective_average({A,-10*A}, 'sphere'))
