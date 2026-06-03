# orthoMTL 0.9.0 (2025-xx-xx)

## Renamed from `orthopen` to `orthoMTL`

This release represents a major refactor and scope expansion of the
`orthopen` package (https://github.com/kevinVervier/orthopen).

## New features

* **Survival analysis mode**: `survival` and `censored.mat` arguments
  in `orthoMTL()` enable censored time-to-event data.
* **Elastic-net sparsity**: `alpha` mixing parameter in [0, 1] blending
  the orthogonality penalty (`alpha = 0`) with L1/Lasso (`alpha = 1`).
* **Survival data utilities**: `create_longitudinal_labels()`,
  `create_indicator_matrix()`, `create_constraint_matrix()`.
* **Cross-validation**: `cv_optimization_foreach()` with parallel
  grid search over lambda, step size, diagonal value, and method.
* **Bootstrap inference**: `bootstrapMTL()` compares real coefficient
  variability against null-model permutations.
* **Evaluation**: `cindex_mtl()` for multi-task concordance index.
* **Visualisation suite**: `plotHeatmap()`, `correlationMap()`,
  `predictionPlot()`, `plot_bootstrapMTL()`.

## Breaking changes from `orthopen`

* Function renamed: `orthopen()` → `orthoMTL()`.
* Return list key renamed: `$W` → `$B`.
* Default `disjoint` changed from `TRUE` to `FALSE`.
* Elastic-net is controlled by a single `alpha` mixing parameter in
  [0, 1] (replacing the `enet`/`lambda1` pair). The penalty is
  `lambda * [(1-alpha)/2 * Omega_K(W)^2 + alpha * ||W||_1]`, ported from
  orthopen v1.1.0. This fixes the previous inconsistent mixing
  (L2 ≈ 0.25*lambda vs L1 = 0.5*lambda1).
* `Iso` package dependency removed; replaced by internal
  `nnmaxheap_C()`.

## Bug fixes

* Variable `T` (masking `base::TRUE`) renamed to `numTasks`.
* Input validation added for `X`, `Y`, `lambda`, `alpha`.
* Warning emitted when `max_iter` is reached without convergence.
* **Logistic gradient corrected** for labels in {-1, +1} (ported from
  orthopen v1.1.0): the gradient previously used the {0, 1} cross-entropy
  form, which is wrong for `y = -1`.
* **Logistic loss made overflow-safe** via the stable softplus identity
  `max(-z,0) + log1p(exp(-|z|))` and exponent clipping in the gradient.

## Infrastructure

* Package scaffolding: LICENSE (GPL-3), `.gitlab-ci.yml`,
  `testthat` runner, `inst/CITATION`.
* Selective NAMESPACE exports (planned for v1.0.0).

---

# orthopen 0.0.1 (2017-01-31)

* Initial release on GitHub.
* Core solver for regression and classification with orthogonal
  column constraints and disjoint support constraints.
* Two vignettes reproducing Figures 2 and 3 from Vervier et al.
  (2014, ECML-PKDD).
