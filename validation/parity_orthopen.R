#!/usr/bin/env Rscript
# =====================================================================
# Parity check: orthoMTL core solver  vs.  orthopen v1.1.0 reference
# ---------------------------------------------------------------------
# orthopen v1.1.0 (github.com/kevinVervier/orthopen) is the INDEPENDENTLY
# VALIDATED solver our REG-01/02/03 fixes were ported FROM. For every case
# the two packages share (everything except orthoMTL's survival masking),
# identical inputs must yield identical coefficient matrices.
#
# The non-survival path is line-for-line identical in both packages
# (zero init, same arithmetic, same stopping rule), so agreement is
# expected to MACHINE PRECISION, not merely "close". Any divergence means
# the port introduced a discrepancy and must be explained.
#
# Run:
#   "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/parity_orthopen.R
# Requires both packages INSTALLED (not just load_all'd).
# =====================================================================

suppressMessages({
  ok_mtl <- requireNamespace("orthoMTL", quietly = TRUE)
  ok_ref <- requireNamespace("orthopen", quietly = TRUE)
})
if (!ok_mtl) stop("orthoMTL is not installed. Run R CMD INSTALL on the repo first.")
if (!ok_ref) {
  # orthopen (github.com/kevinVervier/orthopen) is the reference solver, not a
  # declared dependency. Skip gracefully so a fresh clone / CI does not fail.
  message("SKIP: orthopen (the reference solver) is not installed; ",
          "install it to run the parity check. See validation/README.md.")
  quit(status = 0L)
}

library(orthoMTL)
library(orthopen)

cat(sprintf("orthoMTL %s   vs   orthopen %s (reference)\n\n",
            packageVersion("orthoMTL"), packageVersion("orthopen")))

# Guard: the installed orthoMTL must be the post-migration code (alpha, not enet).
stopifnot("alpha" %in% names(formals(orthoMTL::orthoMTL)))

# ---- shared fixture -------------------------------------------------
set.seed(1)
n <- 80L; p <- 12L; n_tasks <- 3L
X       <- matrix(rnorm(n * p), n, p)
W_true  <- qr.Q(qr(matrix(rnorm(p * n_tasks), p, n_tasks)))
Yreg    <- X %*% W_true + matrix(rnorm(n * n_tasks), n) * 0.1   # regression targets
Ybin    <- ifelse(X %*% W_true > 0, 1, -1)                      # {-1,+1} labels
K       <- matrix(1, n_tasks, n_tasks); diag(K) <- 0.5          # non-trivial task coupling

# common solver controls (identical for both packages so the deterministic
# stopping rule halts both at the same iteration)
ctl <- list(lambda = 1e-2, K = K, step_size = 0.1, tol = 1e-5,
            stop_no_improve = 100, max_iter = 20000)

# ---- helpers --------------------------------------------------------
n_pass <- 0L; n_fail <- 0L
check <- function(label, B_ref, B_new, tol = 1e-8) {
  d <- max(abs(B_ref - B_new))
  pass <- isTRUE(d <= tol)
  if (pass) n_pass <<- n_pass + 1L else n_fail <<- n_fail + 1L
  cat(sprintf("  %-32s max|B_ref - B_new| = %.3e   [%s]\n",
              label, d, if (pass) "PASS" else "FAIL"))
  invisible(pass)
}

ref <- function(Y, ...) do.call(orthopen::orthopen, c(list(X = X, Y = Y), ctl, list(...)))
new <- function(Y, ...) do.call(orthoMTL::orthoMTL, c(list(X = X, Y = Y), ctl, list(...)))

cat("Deterministic cases (zero init, no RNG) -- expect machine precision:\n")

# Case 1 -- regression, orthogonal-columns penalty
check("regression / orthogonal",
      ref(Yreg, disjoint = FALSE)$W,
      new(Yreg, disjoint = FALSE)$B)

# Case 2 -- elastic-net mixing (exercises the alpha L1 term)
check("regression / alpha = 0.5",
      ref(Yreg, disjoint = FALSE, alpha = 0.5)$W,
      new(Yreg, disjoint = FALSE, alpha = 0.5)$B)

check("regression / alpha = 1 (lasso)",
      ref(Yreg, disjoint = FALSE, alpha = 1)$W,
      new(Yreg, disjoint = FALSE, alpha = 1)$B)

# Case 3 -- logistic loss: the REG-01 (gradient convention) + REG-02
#           (overflow-safe softplus) fix path
check("logistic / orthogonal",
      ref(Ybin, logistic = TRUE)$W,
      new(Ybin, logistic = TRUE)$B)

check("logistic / alpha = 0.5",
      ref(Ybin, logistic = TRUE, alpha = 0.5)$W,
      new(Ybin, logistic = TRUE, alpha = 0.5)$B)

cat("\nRNG-seeded case (disjoint support) -- seeds aligned, expect <1e-6:\n")

# Case 4 -- disjoint support uses a random positive init. orthoMTL seeds
# internally (seed = 42); orthopen does not, so we set.seed(42) immediately
# before it to reproduce the identical rnorm() draw, then compare. This also
# exercises proj_disjoint() in both packages.
set.seed(42)
B_ref_dis <- orthopen::orthopen(X, Yreg, lambda = ctl$lambda, K = K,
                                step_size = ctl$step_size, tol = ctl$tol,
                                stop_no_improve = ctl$stop_no_improve,
                                max_iter = ctl$max_iter, disjoint = TRUE)$W
B_new_dis <- orthoMTL::orthoMTL(X, Yreg, lambda = ctl$lambda, K = K,
                                step_size = ctl$step_size, tol = ctl$tol,
                                stop_no_improve = ctl$stop_no_improve,
                                max_iter = ctl$max_iter, disjoint = TRUE, seed = 42)$B
check("regression / disjoint", B_ref_dis, B_new_dis, tol = 1e-6)

# ---- verdict --------------------------------------------------------
cat(sprintf("\n%d passed, %d failed.\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1L)
cat("Core solver matches the validated orthopen v1.1.0 reference.\n")
