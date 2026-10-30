function [node1, node2] = find_max_deg_node(A)
%FIND_MAX_DEG_NODE Find a high-degree adjacent pair of nodes.
%
%   [node1, node2] = FIND_MAX_DEG_NODE(A) finds the node with maximum
%   degree in the graph represented by the adjacency matrix A, then selects
%   among its neighbors the node with maximum degree.
%
%   Input:
%       A       - Adjacency matrix of a graph. Nonzero entries indicate
%                 edges between nodes.
%
%   Outputs:
%       node1   - Index of a node with maximum degree.
%       node2   - Index of the highest-degree neighbor of node1.
%
%   Node degrees are computed as the number of nonzero entries in each row
%   of A. If there are ties, MATLAB's MAX returns the first occurrence.

deg = sum(A ~= 0, 2);        % degree of each node
[~, node1] = max(deg);       % node of maximum degree
N = find(A(node1,:) ~= 0);   % neighbors of node1
[~, ii] = max(deg(N));       % among the neighbors, pick max degree
node2 = N(ii);
end