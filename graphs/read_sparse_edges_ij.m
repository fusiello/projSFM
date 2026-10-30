function [n, ei, ej] = read_sparse_edges_ij(path)
%READ_SPARSE_EDGES_IJ Read an edge list from sparse "(i,j) value" entries.
%
%   [n, ei, ej] = READ_SPARSE_EDGES_IJ(path) reads a text file containing
%   sparse matrix entries of the form
%
%       (i,j) value
%
%   and returns the corresponding undirected edge list.
%
%   Input:
%       path    - Path to the input text file. Lines not matching the
%                 expected sparse-entry format are ignored.
%
%   Outputs:
%       n       - Maximum node index appearing in the parsed entries.
%       ei      - Column vector containing the first endpoint of each edge.
%       ej      - Column vector containing the second endpoint of each edge.
%
%   Only entries whose absolute value is at least 0.5 are retained. Self
%   loops are discarded. Since the graph is treated as undirected, each edge
%   is stored only once, with ei(k) < ej(k). Duplicate edges are removed.
%
%   Node indices in the input file are assumed to be 1-based.


% Parse 1-based "(i,j) 1" lines; return n and undirected unique edges (ei<ej).
pat = '\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*([+-]?\d+(\.\d+)?)';
fid = fopen(path,'r');
if fid < 0, error('Cannot open file: %s', path); end
edges = containers.Map('KeyType','uint64','ValueType','logical');
maxIdx = 0;

while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    tok = regexp(line, pat, 'tokens', 'once');
    if isempty(tok), continue; end
    i = str2double(tok{1});
    j = str2double(tok{2});
    val = str2double(tok{3});
    if abs(val) < 0.5, continue; end
    if i == j, continue; end
    maxIdx = max(maxIdx, max(i,j));
    a = min(i,j); b = max(i,j);
    key = uint64(a) + uint64(maxIdx+1)*uint64(b-1); %#ok<NASGU>
    % Can't use maxIdx in key reliably mid-stream; store as string-like uint64 with large stride:
    key = uint64(a) + uint64(1000000)*uint64(b-1);
    edges(key) = true;
end
fclose(fid);

n = maxIdx;

% extract edges back
ks = cell2mat(keys(edges));
ei = zeros(numel(ks),1);
ej = zeros(numel(ks),1);
for t = 1:numel(ks)
    key = ks(t);
    a = double(mod(key, 1000000));
    b = double(floor(key/1000000) + 1);
    ei(t) = a;
    ej(t) = b;
end

% drop any accidental >n due to stride trick (shouldn't happen if input indices <= 1e6)
ok = (ei>=1 & ei<=n & ej>=1 & ej<=n);
ei = ei(ok); ej = ej(ok);
end