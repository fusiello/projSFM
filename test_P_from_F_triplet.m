%% TEST: triangle reconstruction from fundamental matrices

clear; clc;

rng(0); % reproducibility

tol = 1e-8;

%% 1) Generate random cameras

P1 = random_camera();
P2 = random_camera();
P3 = random_camera();


%% 1) Build a collinear camera triple
% Camera centers lie on the x-axis:
%   C1 = (0,0,0), C2 = (1,0,0), C3 = (2,0,0)

% R1 = eye(3);
% R2 = eye(3);
% R3 = eye(3);
% 
% t1 = [0;0;0];
% t2 = [1;0;0];
% t3 = [2;0;0];
% 
% P1 = [R1, -R1*t1];
% P2 = [R2, -R2*t2];
% P3 = [R3, -R3*t3];

% Ensure they are not degenerate (basic check)
assert(rank(P1)==3 && rank(P2)==3 && rank(P3)==3);

%% 2) Compute fundamental matrices

F12 = fundamental_from_P(P1,P2);
F13 = fundamental_from_P(P1,P3);
F23 = fundamental_from_P(P2,P3);

%% 3) Recover cameras from F

[P1r,P2r,P3r,a,info_rec] = P_from_F_triplet(F12,F13,F23);

% can be used for testing only
if info_rec.degenerate
    disp('Configuration is degenerate.');
end

if ~info_rec.success
    error('Reconstruction failed: %s', info_rec.message);
end

%% 5) Compare fundamental matrices (consistency check)

F12r = fundamental_from_P(P1r,P2r);
F13r = fundamental_from_P(P1r,P3r);
F23r = fundamental_from_P(P2r,P3r);

err12 = projective_distance(F12,F12r);
err13 = projective_distance(F13,F13r);
err23 = projective_distance(F23,F23r);

fprintf('\nFundamental matrix errors:\n');
fprintf('F12 error: %.2e\n', err12);
fprintf('F13 error: %.2e\n', err13);
fprintf('F23 error: %.2e\n', err23);

%% 6) Optional: check projective equivalence

% Solve for H such that P1r ≈ P1 * H
H = relative_projectivity({P1,P2,P3},{P1r,P2r,P3r});

P1_check = P1 * H;
P2_check = P2 * H;
P3_check = P3 * H;

errP1 = projective_distance(P1_check, P1r);
errP2 = projective_distance(P2_check, P2r);
errP3 = projective_distance(P3_check, P3r);

fprintf('\nProjective consistency:\n');
fprintf('P1 error: %.2e\n', errP1);
fprintf('P2 error: %.2e\n', errP2);
fprintf('P3 error: %.2e\n', errP3);

disp('Test completed.');

