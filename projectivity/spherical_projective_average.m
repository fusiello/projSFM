function X = spherical_projective_average(M, weights)
% Frechet spherical mean of a number of unit vectors


weights = weights(:);
s = sum(weights);
if s <= 0
    weights = ones(size(weights)) / numel(weights);
else
    weights = weights / s;
end

sz = size(M{1});
D = prod(sz);
U = zeros(D, numel(M));

% vectorize and normalize
for k = 1:numel(M)
    v = M{k}(:);
    nv = norm(v);
    if nv > 0
        U(:,k) = v / nv;
    end
end

% identify antipodal points
ref = U(:,1);
for k = 2:size(U,2)
    if ref' * U(:,k) < 0
        U(:,k) = -U(:,k);
    end
end

% initialization
c = U * weights;
nc = norm(c);
if nc <= 1e-12
    c = U(:,1);
else
    c = c / nc;
end

for it = 1:1000
    c_prev = c;
    acc = zeros(D,1);

    for k = 1:size(U,2)
        den = sqrt(max(1 - (c_prev' * U(:,k))^2, eps));
        acc = acc + weights(k) * U(:,k) / den;
    end

    na = norm(acc);
    if na <= 1e-12
        break;
    end

    c = acc / na;

    if norm(c - c_prev) < 1e-10
        break;
    end
end

X = reshape(c, sz);
X = X / norm(X, 'fro');
end