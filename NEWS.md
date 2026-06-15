# orthoMTL 0.1.0.9000 (development)

## General-purpose regression & classification support

The solver already supported regression (`logistic = FALSE`,
`survival = FALSE`) and classification (`logistic = TRUE`); these are now
first-class end-to-end with matching evaluation, tuning, and simulation
scaffolding (previously survival-only):

* **Non-survival metrics**: `rmse_mtl()`, `r2_mtl()` (regression) and
  `accuracy_mtl()`, `auc_mtl()` (classification), complementing the
  survival `cindex_mtl()` (unchanged).
* **Mode-aware cross-validation**: `cv_orthoMTL()` gains `logistic` and
  `metric` arguments. The scoring metric now defaults by mode — `cindex`
  for survival, `auc` for logistic, `rmse` otherwise — and selection
  honours each metric's optimisation direction.
* **`predict()` `type` argument**: `"link"` (default, unchanged),
  `"response"` (sigmoid probabilities for logistic fits), and `"class"`
  (predicted {-1, +1} labels for logistic fits).
* **Simulation modes**: `simulate_mtl(mode = ...)` now generates
  `"regression"` and `"classification"` responses in addition to
  `"survival"` (the default).

## Solver

* **Gradient-step schedule exposed** (`schedule` argument to `orthoMTL()`).
  The previously hardcoded `sqrt(i)` decay is now `schedule = "sqrt"` (the
  default, reproducing prior results exactly); `"log"`, `"const"`, and
  `"linear"` are also available. An A/B study
  (`validation/ab_s02_gradient_schedule.R`) found `"log"`/`"const"` reach the
  same optimum ~12--17x faster than `"sqrt"`, while `"linear"` can stall.
  Changing the schedule changes the optimisation path and the exact
  coefficients, so the default is unchanged.

## Behaviour change

* `cv_orthoMTL(survival = FALSE)` previously scored every fold with the
  survival C-index regardless of the data; it now defaults to RMSE for
  plain regression (and AUC when `logistic = TRUE`). Pass `metric =
  "cindex"` to restore the old scoring.

---

# orthoMTL 0.1.0 (2026-06-12)

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
* **Cross-validation**: `cv_orthoMTL()` with parallel grid search over
  lambda, step size, diagonal value, and elastic-net mixing.
* **Bootstrap inference**: `bootstrap_orthoMTL()` compares real
  coefficient variability against null-model permutations.
* **Evaluation**: `cindex_mtl()` for multi-task concordance index.
* **Visualisation suite**: `plot_heatmap()`, `plot_correlation()`,
  `plot_prediction()`, `plot_bootstrap()`.

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

* Package scaffolding: LICENSE (GPL-3), `testthat` runner,
  `inst/CITATION`, GitHub Actions CI.
* Selective NAMESPACE exports.

---

# orthopen 0.0.1 (2017-01-31)

* Initial release on GitHub.
* Core solver for regression and classification with orthogonal
  column constraints and disjoint support constraints.
* Two vignettes reproducing Figures 2 and 3 from Vervier et al.
  (2014, ECML-PKDD).
