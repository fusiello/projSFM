function B = adj2inc(A, oriented)
%ADJ2INC Convert adjacency matrix to sparse incidence matrix.
%
%   B = adj2inc(A, oriented)
%
%   Input:
%     A         - n x n adjacency matrix
%     oriented  - true  -> directed graph
%                 false -> undirected graph
%
%   Output:
%     B         - sparse incidence matrix, size n x m
%
%   Conventions:
%     directed:
%         for edge i -> j, column has -1 at i and +1 at j
%
%     undirected:
%         edges are taken in the order returned by find(triu(A,1))
%         and each column has -1 at the lower-index node
%         and +1 at the higher-index node
%
% same output as  G = graph(A); B = incidence(G)


    if nargin < 2
        oriented = false;
    end

    [n, m] = size(A);
    if n ~= m
        error('A must be square.');
    end

    if oriented
        [src, dst] = find(A);
        keep = (src ~= dst);
        src = src(keep);
        dst = dst(keep);

        ne = numel(src);
        rows = [src; dst];
        cols = [(1:ne)'; (1:ne)'];
        vals = [-ones(ne,1); ones(ne,1)];

    else
        if ~isequal(A, A.')
            error('For oriented = false, A must be symmetric.');
        end

        % Keep original ordering from first version
        [i, j] = find(tril(A,-1));

        % Same ordering, only change signs:
        % -1 on lower index, +1 on higher index
        lo = min(i,j);
        hi = max(i,j);

        ne = numel(i);
        rows = [lo; hi];
        cols = [(1:ne)'; (1:ne)'];
        vals = [-ones(ne,1); ones(ne,1)];
    end

    B = sparse(rows, cols, vals, n, ne);
end