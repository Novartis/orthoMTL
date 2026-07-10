# orthoMTL

> Multi-task learning with orthogonal constraints for time-varying survival analysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/Novartis/orthoMTL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Novartis/orthoMTL/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

Standard survival models such as Cox proportional hazards assume that each
feature's effect is constant over time. **orthoMTL** relaxes this assumption
by recasting survival analysis as a **multi-task learning** problem: binary
classification tasks are defined at multiple time thresholds, and a shared
coefficient matrix is estimated with an orthogonality penalty that encourages
temporally distant tasks to use different features.

The method is based on:

> Vervier, K., Mahé, P., d'Aspremont, A., Veyrieras, J.-B., & Vert, J.-P.
> (2014). On Learning Matrices with Orthogonal Columns or Disjoint Supports.
> *ECML-PKDD 2014*.
> <https://hal.science/hal-00985654>

The survival extension is described in:

> Juric, D. et al. (2026). Genomic Determinants of Response to Alpelisib
> Plus Fulvestrant in the SOLAR-1 Trial. *Annals of Oncology*.
> <https://doi.org/10.1016/j.annonc.2026.04.003>

## Installation

### From CRAN (once released)

```r
install.packages("orthoMTL")
```

### Development version from GitHub

```r
# install.packages("pak")
pak::pak("Novartis/orthoMTL")
```

Or using `remotes`:

```r
# install.packages("remotes")
remotes::install_github("Novartis/orthoMTL")
```

## Quick start

```r
library(orthoMTL)

# Simulate data with known time-varying effects
sim <- simulate_mtl(n = 200, p = 20, n_signals = 5, seed = 42)

# Convert survival times to binary labels at each threshold
Y <- create_longitudinal_labels(sim$SurvTime, sim$Event, sim$thresholds)
W <- create_indicator_matrix(Y)   # censoring mask
K <- create_constraint_matrix(length(sim$thresholds))

# Fit the model
fit <- orthoMTL(sim$X, Y, lambda = 1e-3, K = K,
                survival = TRUE, censored.mat = W)
summary(fit)

# Visualise coefficients across time thresholds
plot_heatmap(fit)
```

### Regression and classification

The same solver handles plain multi-task regression and classification —
survival is just one of three modes. Select the mode on the solver with the
`logistic` / `survival` flags, score with the matching metric, and (for CV)
let `cv_orthoMTL()` pick the metric automatically.

```r
library(orthoMTL)

# --- Multi-task regression ---
sim <- simulate_mtl(n = 200, p = 20, n_signals = 5,
                    mode = "regression", seed = 42)
fit <- orthoMTL(sim$X, sim$Y, lambda = 1e-3)      # logistic = survival = FALSE
pred <- predict(fit, newdata = sim$X)             # type = "link"
rmse_mtl(sim$Y, pred)
r2_mtl(sim$Y, pred)

# --- Multi-task classification (labels in {-1, +1}) ---
sim <- simulate_mtl(n = 200, p = 20, n_signals = 5,
                    mode = "classification", seed = 42)
fit <- orthoMTL(sim$X, sim$Y, lambda = 1e-3, logistic = TRUE)
prob <- predict(fit, newdata = sim$X, type = "response")   # P(Y = +1)
cls  <- predict(fit, newdata = sim$X, type = "class")      # {-1, +1}
auc_mtl(sim$Y, predict(fit, newdata = sim$X))              # AUC on link scores
accuracy_mtl(sim$Y, predict(fit, newdata = sim$X))

# CV picks the metric by mode (cindex / auc / rmse) unless overridden
cv <- cv_orthoMTL(sim$X, sim$Y, logistic = TRUE, survival = FALSE,
                  lambdas = c(1e-3, 1e-2), n_cores = 1)
cv$best
```

For a full walkthrough — including cross-validation, bootstrap inference, and
comparison with Cox — see the package vignette:

```r
vignette("introduction", package = "orthoMTL")
```

## Key functions

| Function | Purpose |
|---|---|
| `orthoMTL()` | Core solver (regression / classification / survival) |
| `cv_orthoMTL()` | Parallel hyperparameter grid search |
| `bootstrap_orthoMTL()` | Bootstrap coefficient stability vs null |
| `simulate_mtl()` | Simulate data with time-varying ground truth |
| `create_longitudinal_labels()` | Survival → binary label matrix |
| `create_indicator_matrix()` | Censoring indicator matrix |
| `create_constraint_matrix()` | Diffusion constraint matrix for tasks |
| `predict()`, `coef()` | Predictions (`type = "link"/"response"/"class"`) and coefficients |
| `cindex_mtl()` | Survival concordance index |
| `rmse_mtl()`, `r2_mtl()` | Regression metrics |
| `accuracy_mtl()`, `auc_mtl()` | Classification metrics |
| `plot_heatmap()`, `plot_correlation()`, `plot_bootstrap()` | Visualisation |

## Citation

If you use orthoMTL in your work, please cite:

```bibtex
@inproceedings{vervier2014,
  author    = {Vervier, Kevin and Mah{\'e}, Pierre and d'Aspremont,
               Alexandre and Veyrieras, Jean-Baptiste and Vert, Jean-Philippe},
  title     = {On Learning Matrices with Orthogonal Columns or Disjoint Supports},
  booktitle = {Machine Learning and Knowledge Discovery in Databases (ECML-PKDD 2014)},
  year      = {2014},
  publisher = {Springer},
  url       = {https://hal.science/hal-00985654}
}
```

## License

GPL-3.
