#!/usr/bin/env Rscript
# =====================================================================
# A/B (S-01): is the loss/penalty normalisation asymmetry a problem?
# ---------------------------------------------------------------------
# KANBAN S-01 (orthoMTL.R:47): the data-fit loss is divided by nrow(X)
# but the penalty terms are not, so "effective regularisation strength
# depends on sample size". This script asks empirically: at a FIXED
# lambda, does the amount of shrinkage drift as n grows?
#
# The objective is
#     A (shipped):   (1/2n) ||XW - Y||^2  +  lambda * Omega_K(W)
# With the default K = I the penalty Omega_K(W) = ||W||_F^2, so this is
# ordinary ridge with an n-AVERAGED loss -- the standard glmnet/cv.glmnet
# convention (lambda is "per-observation").
#
# Dividing the penalty by n as well gives the symmetric variant
#     B (symmetric): (1/2n) ||XW - Y||^2  +  (lambda/n) * Omega_K(W)
# and minimising B at lambda is IDENTICAL to minimising A at lambda/n
# (the 1/n factor is a positive scalar on the whole objective). So we
# realise B by calling the SHIPPED solver with lambda/n -- no solver
# edit, no risk to published numerics.
#
# CRITERION: a sample-size-robust scheme keeps the shrinkage ratio
#   ||B_hat|| / ||B_OLS||  roughly CONSTANT as n grows. A scheme whose
# penalty influence vanishes (or explodes) with n will drift.
#
# Run (against the INSTALLED package):
#   "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/ab_s01_normalisation.R
# =====================================================================

suppressMessages(library(orthoMTL))

set.seed(11)
p        <- 20L
n_tasks  <- 3L
lambda   <- 0.10           # fixed across the whole sweep
noise_sd <- 1.0
n_grid   <- c(40L, 80L, 160L, 320L, 640L, 1280L)

# Fixed ground-truth coefficients and a fixed population covariance so the
# only thing changing across the sweep is the sample size n.
B_true <- matrix(rnorm(p * n_tasks), p, n_tasks)

fro <- function(A) sqrt(sum(A^2))

# Closed-form OLS reference (near-unpenalised) per task, used only to
# normalise the shrinkage ratio. Falls back gracefully if X'X is singular.
ols_fit <- function(X, Y) {
  XtX <- crossprod(X)
  XtX <- XtX + diag(1e-8, ncol(XtX))   # tiny ridge for numerical safety
  solve(XtX, crossprod(X, Y))
}

rows <- vector("list", length(n_grid))

for (i in seq_along(n_grid)) {
  n <- n_grid[i]
  X <- matrix(rnorm(n * p), n, p)
  Y <- X %*% B_true + matrix(rnorm(n * n_tasks, sd = noise_sd), n, n_tasks)

  B_ols <- ols_fit(X, Y)

  # A: shipped normalisation (lambda per-observation)
  fitA <- orthoMTL(X, Y, lambda = lambda, alpha = 0,
                   stop_no_improve = 400L, max_iter = 2e5, seed = 1)
  # B: symmetric normalisation == shipped solver at lambda/n
  fitB <- orthoMTL(X, Y, lambda = lambda / n, alpha = 0,
                   stop_no_improve = 400L, max_iter = 2e5, seed = 1)

  rows[[i]] <- data.frame(
    n          = n,
    convA      = fitA$converged,
    convB      = fitB$converged,
    shrinkA    = fro(fitA$B) / fro(B_ols),   # 1 == no shrinkage
    shrinkB    = fro(fitB$B) / fro(B_ols),
    recovA     = fro(fitA$B - B_true) / fro(B_true),
    recovB     = fro(fitB$B - B_true) / fro(B_true)
  )
}

tab <- do.call(rbind, rows)

cat("A/B S-01: normalisation vs sample size  (p =", p,
    ", tasks =", n_tasks, ", lambda =", lambda, ", K = I -> ridge)\n\n")
cat("  A = shipped  (1/2n)||XW-Y||^2 + lambda*||W||^2   [penalty NOT /n]\n")
cat("  B = symmetric, realised as shipped solver at lambda/n\n\n")

print(format(tab, digits = 3), row.names = FALSE)

# --- A scheme is "good" if its regularisation is (i) ACTIVE -- shrinkage
# meaningfully below 1 -- and (ii) STABLE in n at large n, where the OLS
# denominator is itself stable (the rise at small n is OLS-variance in the
# denominator, not the estimator). We judge stability on the upper half of
# the n grid to avoid that small-sample artefact.
ACTIVE_TOL <- 0.02   # >=2% average shrinkage counts as "actively regularising"
STABLE_TOL <- 0.02   # <=2% range over large-n counts as "stable"

big <- tab$n >= stats::median(tab$n)
regA <- mean(1 - tab$shrinkA); regB <- mean(1 - tab$shrinkB)
rngA <- diff(range(tab$shrinkA[big])); rngB <- diff(range(tab$shrinkB[big]))

cat("\nEffective regularisation summary:\n")
cat(sprintf("      A (shipped)  : avg shrinkage = %.3f  | large-n range = %.3f\n",
            regA, rngA))
cat(sprintf("      B (lambda/n) : avg shrinkage = %.3f  | large-n range = %.3f\n",
            regB, rngB))

A_active <- regA >= ACTIVE_TOL; A_stable <- rngA <= STABLE_TOL
B_active <- regB >= ACTIVE_TOL

cat("\nVERDICT:\n")
if (A_active && A_stable && !B_active) {
  cat("  The shipped per-observation normalisation (A) applies an ACTIVE,",
      "n-STABLE\n  amount of shrinkage across a 32x range of n. The",
      "symmetric variant (B)\n  divides the penalty by n, so its influence",
      "WASHES OUT (shrinkage -> 1, toward\n  OLS): under B you would have to",
      "grow lambda with n to keep regularisation\n  constant -- which is",
      "exactly the sample-size dependence S-01 worries about.\n",
      " It afflicts B, NOT the shipped scheme. A matches the standard glmnet",
      "\n  convention (lambda is per-observation). CONCLUSION: keep the",
      "current\n  normalisation; S-01 is a documented property, not a defect.\n")
} else {
  cat(sprintf("  Mixed/unexpected result (A_active=%s A_stable=%s B_active=%s).",
              A_active, A_stable, B_active),
      "\n  Inspect the table before drawing conclusions about S-01.\n")
}

# Degenerate-result guard (NaNs / non-convergence), not a research verdict.
if (any(!is.finite(unlist(tab[, -1]))) ||
    any(!tab$convA) || any(!tab$convB)) {
  cat("\nWARNING: some fits did not converge or produced non-finite values;",
      "\n  increase max_iter/stop_no_improve before trusting the numbers.\n")
  quit(status = 1L)
}
quit(status = 0L)
