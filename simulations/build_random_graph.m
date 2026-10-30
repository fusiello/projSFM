function [edges, A] = build_random_graph(n, pEdge, addBothDirections, minEdgesPerNode)
%BUILD_RANDOM_GRAPH Generate a random graph with minimum degree and guaranteed connectivity
%
%   [edges, A] = build_random_graph(n, pEdge, addBothDirections, minEdgesPerNode)
%
% INPUT
%   n                  : number of nodes
%   pEdge              : probability of edge creation
%   addBothDirections  : if true, build an undirected graph
%                        if false, build a directed graph
%   minEdgesPerNode    : minimum out-degree per node (or degree if undirected)
%
% OUTPUT
%   edges  : m x 2 edge list
%   A      : n x n sparse adjacency matrix
%
% NOTES
%   - If addBothDirections = true, connectivity is standard undirected connectivity.
%   - If addBothDirections = false, the function guarantees weak connectivity
%     by connecting components while ignoring edge directions.

    if nargin < 3 || isempty(addBothDirections)
        addBothDirections = false;
    end

    if nargin < 4 || isempty(minEdgesPerNode)
        minEdgesPerNode = 0;
    end

    if n < 1
        error('n must be positive.');
    end

    if minEdgesPerNode > n-1
        error('minEdgesPerNode must be <= n-1.');
    end

    if addBothDirections
        % -------------------------------------------------
        % Undirected graph
        % -------------------------------------------------
        M = rand(n) < pEdge;
        M = triu(M, 1);
        M = M | M.';
        M(1:n+1:end) = false;

        deg = full(sum(M, 2));

        for i = 1:n
            while deg(i) < minEdgesPerNode
                candidates = find(~M(i,:));
                candidates(candidates == i) = [];

                if isempty(candidates)
                    break;
                end

                j = candidates(randi(numel(candidates)));

                M(i,j) = true;
                M(j,i) = true;

                deg(i) = deg(i) + 1;
                deg(j) = deg(j) + 1;
            end
        end

        % Guarantee connectivity by connecting components
        G = graph(sparse(M));
        bins = conncomp(G);

        while numel(unique(bins)) > 1
            comps = unique(bins);

            c1 = comps(1);
            c2 = comps(2);

            i = find(bins == c1, 1);
            j = find(bins == c2, 1);

            M(i,j) = true;
            M(j,i) = true;

            G = graph(sparse(M));
            bins = conncomp(G);
        end

    else
        % -------------------------------------------------
        % Directed graph
        % -------------------------------------------------
        M = rand(n) < pEdge;
        M(1:n+1:end) = false;

        degOut = full(sum(M, 2));

        for i = 1:n
            while degOut(i) < minEdgesPerNode
                candidates = find(~M(i,:));
                candidates(candidates == i) = [];

                if isempty(candidates)
                    break;
                end

                j = candidates(randi(numel(candidates)));

                M(i,j) = true;
                degOut(i) = degOut(i) + 1;
            end
        end

        % Guarantee weak connectivity by connecting weak components
        Gweak = graph(sparse(M | M.'));
        bins = conncomp(Gweak);

        while numel(unique(bins)) > 1
            comps = unique(bins);

            c1 = comps(1);
            c2 = comps(2);

            i = find(bins == c1, 1);
            j = find(bins == c2, 1);

            % Add one directed edge between components
            M(i,j) = true;

            Gweak = graph(sparse(M | M.'));
            bins = conncomp(Gweak);
        end
    end

    A = spones(sparse(M));
    [ii, jj] = find(A);
    edges = [ii, jj];
end