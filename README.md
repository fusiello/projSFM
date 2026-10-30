# projSFM — MATLAB tools for projective structure from motion

Helper functions for projective structure from motion: distances, averages and
relative projectivities between projective elements, camera recovery from
fundamental matrices, and finite-solvability tests on viewing graphs.

## Requirements

- MATLAB R2017b or later (uses `graph`, `vecnorm`, string literals and local
  functions in scripts). No additional toolboxes are needed.
- The code is not Octave-compatible as a whole (Octave lacks the `graph` class
  and some `eigs`/`svds` options).

## Getting started

```matlab
cd projSFM
setup_projSFM        % adds the sub-folders to the MATLAB path
test_solvability     % runs the solvability tests on the graphs in data/
```

## Folder layout

| Folder | Contents |
|---|---|
| `projectivity/` | Distances, averages and relative projectivities between projective elements |
| `basicCV/` | Basic multi-view geometry (F from P, P from a triplet of F, small helpers) |
| `solvability/` | Construction of solvability matrices and numerical nullity estimation — see [`solvability/README.md`](solvability/README.md) |
| `graphs/` | Adjacency / incidence / edge-list utilities |
| `simulations/` | Random cameras, graphs and perturbations for synthetic experiments |
| `data/` | Example 90-node viewing graphs (adjacency matrices) used by `test_solvability` |

### `projectivity/`

| File | Description |
|---|---|
| `projective_distance.m` | Distance between two projective elements (`'sphere'` or `'op'` metric) |
| `projective_average.m` | Average of projective elements (calls `spherical_projective_average`) |
| `spherical_projective_average.m` | Average of projective elements on the sphere |
| `relative_projectivity.m` | 4×4 projectivity between two sets of 3×4 or 4×4 matrices (least squares) |
| `relative_projectivity_l1.m` | Robust L1 (DPCP-IRLS, Huber-smoothed) version of `relative_projectivity` |

### `basicCV/`

| File | Description |
|---|---|
| `fundamental_from_P.m` | F from two cameras (`'minors'` or `'pseudo'`) |
| `P_from_F_triplet.m` | Recover three consistent cameras from the three Fs of a triangle (à la Colombo) |
| `inv3x3_fast.m` | Closed-form inverse of a 3×3 matrix |
| `fund.m`, `skew.m`, `commutation.m` | Minimal stand-ins for the Computer Vision Toolkit functions of the same name, so the package is self-contained |

### `graphs/`

| File | Description |
|---|---|
| `adj2inc.m` | Adjacency → incidence matrix |
| `inc2edgelist.m` | Incidence matrix → edge list |
| `find_max_deg_node.m` | Highest-degree node and its highest-degree neighbour |
| `read_sparse_edges_ij.m` | Read an undirected edge list from a text dump of sparse `(i,j) value` entries |

### `simulations/`

| File | Description |
|---|---|
| `build_random_graph.m` | Random connected graph |
| `random_camera.m` | Random full-rank camera |
| `random_cams_funds.m` | Random cameras (and fundamental matrices) on a graph |
| `funds_from_cams.m` | Fundamental matrices on the edges of a graph |
| `perturb_by_angle_on_sphere.m` | Angular perturbation of a projective element |

## Test scripts

| Script | What it checks |
|---|---|
| `test_distance_and_average.m` | Scale/sign invariance of distance and average |
| `test_relative_projectivity.m` | Recovery of a known projectivity |
| `test_P_from_F_triplet.m` | Camera recovery from a triangle of fundamental matrices |
| `test_solvability.m` | `solvability_all` on the graphs in `data/`, all methods × all modes |

## License

Released under the MIT License (see [`LICENSE`](LICENSE)): you may use, modify
and redistribute the code, provided the copyright notice is kept.
`solvability/ParallelRigidityTest.m` is adapted from code by Onur Ozyesil and
remains subject to its author's terms.

## Citation

If you use this code in academic work, please cite it (see
[`CITATION.cff`](CITATION.cff)):

```
A. Fusiello. projSFM: MATLAB tools for projective structure from motion.
University of Udine, 2026.
```
