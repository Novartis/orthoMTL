#!/usr/bin/env Rscript
# =====================================================================
# A/B (S-04): are the convergence defaults enough at high task count?
# ---------------------------------------------------------------------
# KANBAN S-04 (orthoMTL.R:55): the default stop_no_improve = 100 (and,
# secondarily, max_iter) may be too tight for fine-grained survival
# problems. Prior evidence: the masked survival solver hit the iteration
# cap at k = 16 threshold-tasks while converging in < 1k iters at k = 4.
#
# This script holds the DATA and lambda fixed and varies only the
# convergence budget, across an increasing number of survival tasks k:
#     A (shipped patience): stop_no_improve = 100
#     B (generous budget) : stop_no_improve = 1000
# Both use the same (capped) max_iter so a fit that "wants" more steps is
# visible as hitting the cap. No solver edit -- both are plain arguments.
#
# We report, per k: did each fit converge, how many iterations it used,
# the final objective, and the held-out multi-task C-index. If B reaches
# a materially lower objective / higher C-index than A at large k, the
# default patience is too tight there.
#
# Run (against the INSTALLED package):
#   "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/ab_s04_convergence.R
# =====================================================================

suppressMessages(library(orthoMTL))

# A deliberately HARD regime: p approaches n (under-determined) and a small
# lambda, where convergence is slowest, so any inadequacy of the defaults
# would show. One population, re-labelled at k thresholds per row.
set.seed(23)
n_tr   <- 120L
n_te   <- 200L
p      <- 40L
max_it <- 50000L          # shared cap; "conv == FALSE" means it hit this cap
lam_grid <- c(1e-2, 1e-4) # 1e-4 is the slow-to-converge end
k_grid   <- c(4L, 8L, 12L, 16L)

sim <- simulate_mtl(n = n_tr + n_te, p = p, n_signals = 6,
                    mode = "survival", seed = 23)
tr <- seq_len(n_tr); te <- n_tr + seq_len(n_te)

fit_one <- function(X, Y, W, lambda, patience) {
  orthoMTL(X, Y, lambda = lambda, survival = TRUE, censored.mat = W,
           stop_no_improve = patience, max_iter = max_it, seed = 1, verbose = 0)
}

rows <- list()
for (lambda in lam_grid) {
  for (k in k_grid) {
    # k evenly spaced interior quantiles of the observed times -> k tasks.
    thr <- as.numeric(stats::quantile(sim$SurvTime, probs = seq_len(k) / (k + 1)))
    thr <- sort(unique(round(thr, 6)))
    Yall <- create_longitudinal_labels(sim$SurvTime, sim$Event, thr)
    Wall <- create_indicator_matrix(Yall)

    Xtr <- sim$X[tr, ]; Ytr <- Yall[tr, ]; Wtr <- Wall[tr, ]
    Xte <- sim$X[te, ]; Yte <- Yall[te, ]

    fitA <- fit_one(Xtr, Ytr, Wtr, lambda, patience = 100L)
    fitB <- fit_one(Xtr, Ytr, Wtr, lambda, patience = 1000L)

    predA <- t(apply(Xte %*% fitA$B, 1, nnmaxheap_C))
    predB <- t(apply(Xte %*% fitB$B, 1, nnmaxheap_C))

    rows[[length(rows) + 1L]] <- data.frame(
      lambda  = lambda,        k       = length(thr),
      A_iter  = fitA$imax,     A_obj   = fitA$obj,  A_cidx = cindex_mtl(Yte, predA),
      B_iter  = fitB$imax,     B_obj   = fitB$obj,  B_cidx = cindex_mtl(Yte, predB),
      hitcap  = (!fitA$converged) || (!fitB$converged),
      objdrop = (fitA$obj - fitB$obj) / abs(fitA$obj)  # >0 => B (patient) lower
    )
  }
}
tab <- do.call(rbind, rows)

cat("A/B S-04: convergence budget vs task count  (n_tr =", n_tr,
    ", p =", p, ", max_iter =", max_it, ")\n\n")
cat("  Two levers tested independently:\n")
cat("   - PATIENCE: A stop_no_improve=100  vs  B stop_no_improve=1000\n")
cat("   - MAX_ITER: shared cap; A_iter/B_iter << cap means the cap never binds\n\n")

print(format(tab, digits = 4), row.names = FALSE)

# Lever 1: does patience change the SOLUTION? (objdrop ~ 0 => no)
max_objdrop <- max(abs(tab$objdrop))
# Lever 2: does the iteration cap ever bind?
cap_binds   <- any(tab$hitcap)
max_iter_used <- max(tab$A_iter, tab$B_iter)

cat(sprintf("\nPatience lever:  max |objective change A->B| = %.2e\n", max_objdrop))
cat(sprintf("Max_iter lever:  most iterations used = %d of %d cap  (cap binds: %s)\n",
            max_iter_used, max_it, cap_binds))
cat("Iteration scaling: iters grow with k and 1/lambda (see A_iter column).\n")

cat("\nVERDICT:\n")
MAT_OBJ <- 0.01
if (max_objdrop <= MAT_OBJ && !cap_binds) {
  cat("  Neither lever is the bottleneck in this hard regime:\n",
      " - extra patience leaves the objective unchanged (the no-improve rule",
      "is NOT\n    stopping early -- it finds the same optimum), and\n",
      " - the iteration count, though it grows with k and 1/lambda, stays far",
      "below\n    the cap (and the shipped default max_iter = 1e6 is larger",
      "still).\n",
      " The current defaults are adequate for the post-REG-fix solver. The",
      "earlier\n  'k=16 hits 20k cap' evidence likely predates those fixes or",
      "used a low\n  max_iter. Recommend: keep stop_no_improve = 100; if",
      "anything, document that\n  max_iter (not patience) is the lever for very",
      "fine-grained / low-lambda fits.\n")
} else if (max_objdrop > MAT_OBJ) {
  cat(sprintf("  Patience IS a lever: more patience lowers the objective by up",
      "to %.1f%%.\n  Recommend raising stop_no_improve (or scaling it with",
      "k). S-04 confirmed.\n", 100 * max_objdrop))
} else {
  cat("  The iteration cap binds in this regime: max_iter is the lever, not",
      "patience.\n  Recommend documenting/raising max_iter for fine-grained",
      "survival. S-04 (partial).\n")
}

if (any(!is.finite(unlist(tab)))) {
  cat("\nWARNING: non-finite values in the table; inspect before trusting.\n")
  quit(status = 1L)
}
quit(status = 0L)
