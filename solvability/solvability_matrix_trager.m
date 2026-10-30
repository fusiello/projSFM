function [S]=solvability_matrix_trager(G)
% Original version from the paper:
% VIEWING GRAPH SOLVABILITY IN PRACTICE
% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023
%
% Build the reduced finite-solvability matrix and augment it with
% gauge-fixing equations so that the target is full column rank.

n = numnodes(G);
m = numedges(G);       % number of edges in input graph
deg = degree(G);       % degree of nodes

C = rand(4,numnodes(G))*100; % TO CHECK: which is the best range?

%% build the REDUCED finite solvability matrix

disp('Building the reduced solvability matrix...')

nEq = sum(11*(deg-1)); % number of equations

A = spalloc(nEq,16*m,nEq*16*2); % initialization (#rows, #columns, #nonzero)

k=1; % index current equation
for i=1:n % for all the nodes ...
    % let us consider node i in the input graph

    % build local equations (depend on camera centers)

    Ai=finite_equations_independent(C(:,i)); % 20 equations (redundant)


    N = neighbors(G,i);  % find nieghbors of node i
    assert(length(N)==deg(i)) % sanity check

    tau=findedge(G,i,N(1)); % first neighbor is the root
    for j=2:deg(i)
        nu=findedge(G,i,N(j)); % pick another edge in the tree

        % update the global matrix

        A(11*k-10:11*k,16*tau-15:16*tau)=Ai;
        A(11*k-10:11*k,16*nu-15:16*nu)=-Ai;


        k=k+1;
    end

end


%% Fix projective & scale ambiguity by adding rows

disp('Fixing the trivial d.o.f.')

% fix first transformation to be something
B=sparse(16,16*m);
B(1:16,1:16)=speye(16);

% faster version
i     = repmat(1:(m-1),16,1);
D     = sparse(i(:), 17:m*16, ones(1,16*(m-1)));

S=[A;B;D];

end



%% VIEWING GRAPH SOLVABILITY IN PRACTICE
%% Federica Arrigoni, Tomas Pajdla, Andrea Fusiello. ICCV 2023

function A=finite_equations_independent(c)
% INPUT: camera centre
% output: 11x16 matrix associated with the camera centre, which is the
% building block of the REDUCED finite solvability equations


% Legacy Federica's version
% A=finite_equations(c);
% % select only independent equations (11 out of 20)
% A=A([1:5 7:9 13:15],:);


C1=c(1);
C2=c(2);
C3=c(3);
C4=c(4);


i = [11
    7
    6
    8
    4
    10
    3
    9
    5
    9
    2
    1
    10
    9
    11
    5
    8
    10
    2
    4
    11
    1
    7
    3
    6
    11
    9
    10];


j = [ 1
    2
    3
    4
    5
    6
    7
    7
    8
    8
    9
    10
    10
    11
    11
    12
    12
    12
    13
    13
    13
    14
    14
    15
    15
    15
    16
    16];


v =[-C4
    -C4
    -C4
    -C3
    -C4
    -C3
    -C4
    -C3
    -C3
    C4
    -C4
    -C4
    C2
    C2
    C4
    C2
    C1
    -C4
    C3
    C2
    C1
    C3
    C1
    C2
    C1
    -C3
    -C2
    C3];

A = sparse(i,j,v,11,16);

end