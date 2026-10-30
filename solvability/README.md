# Solvability Matrix Construction and Nullity Estimation

This folder contains MATLAB code for constructing solvability matrices associated with several graph-based finite-solvability tests, and for estimating their numerical nullity from small singular values.

The main entry point is:

```matlab
summary = solvability_all(A, kSmall, tol, method, mode);
```

where `A` is the adjacency matrix of the input graph. Note that `method` comes **before** `mode`.

---

## Files

### `solvability_all.m`

Main wrapper that builds one of the available solvability matrices and estimates its numerical nullity.

Supported construction methods:

```matlab
method = 'trager' % ICCV23 reduced Trager-Osserman-Ponce construction
method = 'pami'   % legacy ECCV/PAMI skew-symmetry construction
method = 'skew'   % optimized skew-symmetry construction
method = 'fund'   % optimized fundamental-matrix construction
```

Typical usage:

```matlab
summary = solvability_all(A, 6, 5e-8, 'skew', 'eigs_sm');
```

The returned `summary` structure contains fields such as:

```matlab
summary.nullity_est
summary.s_small
summary.smax
summary.J
summary.n
summary.m_undirected
summary.eq
summary.vars
summary.thr
summary.kUse
```

The graph is considered numerically solvable when

```matlab
summary.nullity_est < 1
```

---

### `estimate_nullity.m`

Estimates the numerical nullity of a matrix from its smallest singular values.

Typical usage:

```matlab
[nullity_est, s_small, stats] = estimate_nullity(J, k, tol, mode);
```

Supported modes include:

```matlab
mode = 'svd'      % full dense SVD
mode = 'svds'     % sparse smallest singular values
mode = 'eigs_sm'  % eigenvalue estimate from J'*J, smallest magnitude
mode = 'eigs'     % eigenvalue estimate from J'*J, shift-invert
mode = 'auto'     % try eigs_sm, then fallback to shift-invert
mode = 'gt_nulls' % null-space dimension via gt_nulls (Gotsman-Toledo)
```

The estimated nullity is

```matlab
nullity_est = sum(s_small < tol);
```

The input matrix is assumed to be row-normalized by the caller when required.

---

### `gt_nulls.m`

Null space of a (possibly sparse) rectangular matrix, after Gotsman and Toledo,
"On the computation of null spaces of sparse rectangular matrices".

```matlab
[N, flag, bound] = gt_nulls(A, k, tol);
```

`N` has orthonormal null vectors as columns; `flag == 0` means `bound` is the
exact null-space dimension, `flag == 1` means it is an upper bound. Used by
`estimate_nullity` in `'gt_nulls'` mode and by `ParallelRigidityTest`.

---

### `ParallelRigidityTest.m`

Parallel-rigidity test of a graph in R^d (after O. Ozyesil), optionally
returning the largest maximally parallel-rigid component.

```matlab
[isRigid, Acomp, testInfo] = ParallelRigidityTest(A, d, tolNull, mode);
```

---

### `solvability_matrix_trager_a.m`, `solvability_matrix_trager.m`

Reduced Trager–Osserman–Ponce construction (ICCV23); `_a` is the optimized
version used by `solvability_all(..., 'trager')`. The input is a MATLAB
`graph` object.

---

### `solvability_matrix_fund.m`

Original (unoptimized) fundamental-matrix construction, kept for reference;
see `solvability_matrix_fund_a.m`.

---

### `solvability_matrix_skew.m`

Legacy skew-symmetry solvability-matrix construction.

This is the original ECCV/PAMI-style implementation. It accepts a MATLAB `graph` object and builds the finite solvability matrix from all graph edges.

Typical usage:

```matlab
[A,m,Cams,Funds] = solvability_matrix_skew(G, true, true, true, node1, node2);
```

The output matrix includes:

- edge-based skew-symmetry constraints,
- optional projective gauge constraints,
- optional scale constraints,
- optional full-rank camera constraints.

This version samples cameras and fundamental matrices internally.

---

### `solvability_matrix_skew_a.m`

Optimized skew-symmetry solvability-matrix construction.

Typical usage:

```matlab
S = solvability_matrix_skew_a(A);
```

The input `A` is an adjacency matrix. The function samples cameras and fundamental matrices internally, then builds a sparse matrix using a reduced block pattern to avoid inserting known zero entries.


---

### `solvability_matrix_fund_a.m`

Optimized fundamental-matrix solvability-matrix construction.

Typical usage:

```matlab
Js = solvability_matrix_fund_a(A);
```

The input `A` is an adjacency matrix. The function samples cameras internally and assembles the normalized fundamental-matrix Jacobian.

Each graph edge contributes 8 equations. Gauge columns are removed during assembly, so for a graph with `n` vertices and `m` undirected edges,

```matlab
size(Js,1) = 8*m
size(Js,2) = 11*n - 16
```

The returned matrix is not row-normalized.

---

#
## Basic Usage

```matlab
S = load('data/A140.mat');   % 90-node example graph (run from the projSFM root)
f = fieldnames(S);
A = S.(f{1});

kSmall = 6;
tol    = 5e-8;
method = 'skew';
mode   = 'eigs_sm';

summary = solvability_all(A, kSmall, tol, method, mode);
disp(summary.nullity_est);
```


---

## Input Graphs

Most functions expect the graph as an adjacency matrix `A`.

The matrix should be:

- square,
- symmetric,
- sparse or dense,
- interpreted by its nonzero pattern.

For example:

```matlab
A = sparse(A ~= 0);
A = spones(A | A.');
A = A - diag(diag(A));
```

For the legacy `solvability_matrix_skew.m`, the input is instead a MATLAB `graph` object:

```matlab
G = graph(A);
```

---

## Row Normalization

`solvability_all.m` row-normalizes the constructed solvability matrix before calling `estimate_nullity`.

The normalization is performed with a guard against zero rows:

```matlab
nr = vecnorm(S,2,2);
nr(nr == 0) = 1;
S = S ./ nr;
```

Individual matrix-construction functions may return unnormalized matrices.

---

## Numerical Nullity

The numerical nullity is estimated by counting how many of the computed small singular values are below a threshold:

```matlab
nullity_est = sum(s_small < tol);
```

The value of `kSmall` should be chosen large enough to capture the expected nullity. If the true nullity is larger than `kSmall`, the returned estimate is truncated by `kSmall`.

For robust testing, use for example:

```matlab
kSmall = 6;
tol = 5e-8;
```

---

## Suggested Test

The script `test_solvability.m` in the parent folder runs `solvability_all`
on the graphs in `data/` for every method and mode, and compares the
estimated nullity against the expected result:

```matlab
G90-139 (data/A139.mat) : summary.nullity_est > 0    % not solvable
G90-140 (data/A140.mat) : summary.nullity_est == 0   % solvable
G90-170 (data/A170.mat) : summary.nullity_est == 0   % solvable
```

---

## Notes

- `mode = 'svd'` can be slow or memory intensive because it forms a dense matrix.
- `mode = 'eigs_sm'` and `mode = 'eigs'` are usually more appropriate for large sparse matrices.
- The `pami` method is kept mainly for comparison with the legacy implementation.
- The `skew`, `fund`, and `4dof` methods are the current optimized construction variants.
- `pami`, `trager` and `4dof` use MATLAB `graph` objects internally.
