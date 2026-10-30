function [Js]=solvability_matrix_fund(A)
% Original Irena's version

disp('Building the solvability matrix (fund) ...')

% convert adjacency
B = adj2inc(A, false);
E = inc2edgelist(B);
[n,m] =  size(B);

[cams,~]=random_cams_funds(E,n,m);

rows=[];columns=[];values=[];
for k=1:m % for all the edges ...

    % let us consider edge k in the input graph
    i=E(k,1); % left node
    j=E(k,2); % right node

    P1=cams(:,:,i); % left camera
    P2=cams(:,:,j); % right camera
    %F=funds(:,:,k); % fundamental matrix

    V1=P1(:,1:3);
    m1=P1(:,4);
    V2=P2(:,1:3);
    m2=P2(:,4);
    [JF_P1_c, JF_P2_c] = fundamental_jacobian_C(V1,m1,V2,m2);
    J1=JF_P1_c(:,1:11);
    J2=JF_P2_c(:,1:11);


    rows=[rows repelem(8*(k-1)+1:8*(k-1)+8,11) repelem(8*(k-1)+1:8*(k-1)+8,11)];
    columns=[columns repmat(11*(i-1)+1:11*(i-1)+11,1,8) repmat(11*(j-1)+1:11*(j-1)+11,1,8)];
    J1n=J1.';J2n=J2.';
    values=[values J1n(:)' J2n(:)'];

end

Js=sparse(rows,columns,values,8*m,11*n);
Js(:,[1:12 15 18 21])=[];

end

function [JF_P1, JF_P2] = fundamental_jacobian_C(V1, m1, V2, m2)
    
    P1 = [V1, m1]; 
    P2 = [V2, m2];
    invV1 = inv(V1);
    invV2 = inv(V2);
    c1 = -invV1*m1;
    c2 = -invV2*m2;
    
    e21 = m2-V2*invV1*m1;
    e21_skew = [ 0    -e21(3)  e21(2); 
                 e21(3)  0    -e21(1); 
                -e21(2)  e21(1)  0 ];
    P1_plus = P1' * inv(P1 * P1');

    % skew lifting matrix
    S = [ 0  0  0 ;
          0  0  1 ;
          0 -1  0 ;
          0  0 -1 ;
          0  0  0 ;
          1  0  0 ;
          0  1  0 ;
         -1  0  0 ;
          0  0  0 ];


    % ---- P1 JACOBIAN -----------------------

    I3 = eye(3);
    I4 = eye(4);
    De21_V1 = kron((invV1 * m1).', V2 * invV1);
    De21_m1 = -V2 * invV1;
    De21_P1 = [De21_V1, De21_m1];
    J_P1_A = kron((P2 * P1_plus).', I3) * S * De21_P1;

    % K34 = zeros(12,12);
    % idx = reshape(1:12, 3, 4);
    % idx_T = idx';
    % for i = 1:12
    %     K34(i, idx_T(i)) = 1;
    % end
    K34=commutation(3,4);
    
    DP1p = -kron(P1_plus.', P1_plus)+...
           (kron((I3 - P1*P1_plus), P1_plus*P1_plus.') + ...
            kron((P1_plus.'*P1_plus), (I4 - P1_plus*P1))) * K34;
    J_P1_B = kron(I3, e21_skew * P2) * DP1p;

    JF_P1_raw = J_P1_A + J_P1_B;


    % --- P2 JACOBIAN ---

    De21_V2 = -kron((invV1 * m1).', eye(3));
    De21_m2 = eye(3);
    De21_P2 = [De21_V2, De21_m2];
    J_P2_A = kron((P2 * P1_plus).', eye(3)) * S * De21_P2;

    J_P2_B = kron(P1_plus.', e21_skew);

    JF_P2_raw = J_P2_A + J_P2_B;

     % --- normalize F ---
    F = e21_skew*P2*P1_plus;
    F33 = F(3,3);
    f = F(:);
    e33 = zeros(9,1); e33(9) = 1;

    % Derivative of normalization
    N = (1/F33)*eye(9) - (1/F33^2) * f * e33.';

    JF_P1n = N * JF_P1_raw;
    JF_P2n = N * JF_P2_raw;

    % --- affine chart: drop 9th row ---
    SF = eye(9); SF(9,:) = [];
    JF_P1 = SF * JF_P1n;
    JF_P2 = SF * JF_P2n;
   
end