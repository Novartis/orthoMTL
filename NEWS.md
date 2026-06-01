# orthoMTL 0.9.0 (2025-xx-xx)

## Renamed from `orthopen` to `orthoMTL`

This release represents a major refactor and scope expansion of the
`orthopen` package (https://github.com/kevinVervier/orthopen).

## New features

* **Survival analysis mode**: `survival` and `censored.mat` arguments
  in `orthoMTL()` enable censored time-to-event data.
* **Elastic-net sparsity**: independent `lambda1` parameter for L1
  regularisation on top of the orthogonal penalty.
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
* Elastic-net `lambda1` is now an explicit scalar parameter
  (was internally derived as a vector coupled to `diag(K)`).
* `Iso` package dependency removed; replaced by internal
  `nnmaxheap_C()`.

## Bug fixes

* Variable `T` (masking `base::TRUE`) renamed to `numTasks`.
* Input validation added for `X`, `Y`, `lambda`, `lambda1`.
* Warning emitted when `max_iter` is reached without convergence.

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
