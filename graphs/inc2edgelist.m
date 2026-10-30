function E = inc2edgelist(B)
%INC2EDGELIST Convert a signed incidence matrix to an edge list.
%
%   E = INC2EDGELIST(B) converts the signed incidence matrix B into an
%   edge list E. Each column of B is assumed to represent one directed edge,
%   with exactly one entry equal to -1 and exactly one entry equal to +1.
%
%   Input:
%       B   - Signed incidence matrix. Each column represents one edge and
%             must contain exactly one -1 and one +1.
%
%   Output:
%       E   - m-by-2 edge list, where m = size(B,2). Each row E(k,:) =
%             [src dst] gives the source and destination nodes of edge k.
%
%   The source node is taken as the row containing -1, and the destination
%   node is taken as the row containing +1. The function checks that each
%   column satisfies this signed-incidence structure.

[i,j,v] = find(B);
neg = v < 0;
pos = v > 0;
src = i(neg);
dst = i(pos);
col_neg = j(neg);
col_pos = j(pos);
if numel(src) ~= size(B,2) || numel(dst) ~= size(B,2) || ~isequal(col_neg, col_pos)
    error('B must have exactly one -1 and one +1 per column.');
end
E = [src, dst];
end