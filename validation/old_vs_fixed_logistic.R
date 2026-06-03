#!/usr/bin/env Rscript
# =====================================================================
# Old-vs-fixed A/B: the REG-01 logistic gradient bug, made visible
# ---------------------------------------------------------------------
# The parity check (validation/parity_orthopen.R) shows orthoMTL matches
# the FIXED orthopen v1.1.0 bit-for-bit. This script shows the other half:
# that the fix actually MATTERS, by comparing against the BUGGY pre-v1.0.1
# ancestor (orthopen/new/legacy_orthopen.R) that orthoMTL was forked from.
#
# The bug (REG-01): the legacy logistic LOSS is the correct y in {-1,+1}
# softplus  log(1 + exp(-y*Xw))  (legacy line 99), but its GRADIENT is the
# y in {0,1} cross-entropy gradient  sigma(Xw) - y  (legacy lines 139/155/
# 157). Loss and gradient use different label conventions -> the solver
# descends the wrong direction for y = -1.
#
# RIGOROUS CRITERION (no ground truth needed): the logistic objective
#   sum_k  log(1 + exp(-Y * Xw))  +  penalty(W)
# is invariant under (Y, W) -> (-Y, -W). Hence ANY correct minimiser must
# satisfy  B(-Y) = -B(Y)  (label-flip antisymmetry). The fixed solver
# satisfies it; the buggy one cannot, because its gradient breaks the
# symmetry. We also confirm the fix changes the fit and improves held-out
# logistic loss.
#
# Run:
#   "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/old_vs_fixed_logistic.R
# Requires orthoMTL INSTALLED; reads the legacy source directly (not installed).
# =====================================================================

suppressMessages(library(orthoMTL))

# The buggy pre-v1.0.1 ancestor is read directly from the sibling orthopen
# repo (not installed, not a dependency). Allow override via env var so the
# path is not hard-wired to one machine; skip gracefully if unavailable.
legacy_path <- Sys.getenv("ORTHOPEN_LEGACY",
                          "C:/Users/Shadow/Bricoles/orthopen/new/legacy_orthopen.R")
if (!file.exists(legacy_path)) {
  message("SKIP: legacy (buggy) orthopen source not found at '", legacy_path,
          "'. Set ORTHOPEN_LEGACY to its path to run this A/B. ",
          "See validation/README.md.")
  quit(status = 0L)
}
leg_env <- new.env()
suppressMessages(sys.source(legacy_path, envir = leg_env))
orthopen_buggy <- get("orthopen", envir = leg_env)   # pre-v1.0.1, REG-01 present

# ---- logistic fixture (noisy, not perfectly separable) --------------
set.seed(7)
n <- 160L; p <- 10L; n_tasks <- 2L
X     <- matrix(rnorm(n * p), n, p)
W_t   <- matrix(rnorm(p * n_tasks), p, n_tasks)
prob  <- 1 / (1 + exp(-(X %*% W_t)))
Y     <- ifelse(matrix(runif(n * n_tasks), n, n_tasks) < prob, 1, -1)
K     <- diag(1, n_tasks)
lambda <- 1e-2

tr  <- 1:120; te <- 121:160
Xtr <- X[tr, ]; Ytr <- Y[tr, ]
Xte <- X[te, ]; Yte <- Y[te, ]

# stable mean logistic loss  mean log(1 + exp(-y*Xw))
logloss <- function(B, Xm, Ym) { z <- Ym * (Xm %*% B); mean(pmax(-z, 0) + log1p(exp(-abs(z)))) }
maxabs  <- function(A) max(abs(A))

# ---- fits: fixed (orthoMTL) and buggy (legacy), on Y and on -Y -------
fx_pos <- orthoMTL(Xtr,  Ytr, lambda = lambda, K = K, logistic = TRUE, disjoint = FALSE)$B
fx_neg <- orthoMTL(Xtr, -Ytr, lambda = lambda, K = K, logistic = TRUE, disjoint = FALSE)$B
bg_pos <- orthopen_buggy(Xtr,  Ytr, lambda = lambda, K = K, logistic = TRUE, disjoint = FALSE)$W
bg_neg <- orthopen_buggy(Xtr, -Ytr, lambda = lambda, K = K, logistic = TRUE, disjoint = FALSE)$W

# ---- (1) label-flip antisymmetry  B(-Y) == -B(Y) --------------------
asym_fixed <- maxabs(fx_pos + fx_neg)   # 0 if the symmetry holds
asym_buggy <- maxabs(bg_pos + bg_neg)

cat("Old-vs-fixed logistic A/B  (orthoMTL  vs  pre-v1.0.1 legacy orthopen)\n\n")
cat("(1) Label-flip antisymmetry  max|B(Y) + B(-Y)|   (0 = correct):\n")
cat(sprintf("      fixed (orthoMTL) : %.3e   %s\n", asym_fixed,
            if (asym_fixed < 1e-6) "[holds]" else "[VIOLATED]"))
cat(sprintf("      buggy (legacy)   : %.3e   %s\n", asym_buggy,
            if (asym_buggy < 1e-6) "[holds]" else "[VIOLATED]"))

# ---- (2) the fix changes the fit ------------------------------------
cat(sprintf("\n(2) Fix changes the solution   max|B_fixed - B_buggy| = %.3e\n",
            maxabs(fx_pos - bg_pos)))

# ---- (3) held-out logistic loss (lower = better) --------------------
ll_fixed <- logloss(fx_pos, Xte, Yte)
ll_buggy <- logloss(bg_pos, Xte, Yte)
cat(sprintf("\n(3) Held-out mean logistic loss (lower = better):\n"))
cat(sprintf("      fixed (orthoMTL) : %.4f\n", ll_fixed))
cat(sprintf("      buggy (legacy)   : %.4f   (%+.1f%% vs fixed)\n",
            ll_buggy, 100 * (ll_buggy - ll_fixed) / ll_fixed))

# ---- verdict --------------------------------------------------------
ok <- isTRUE(asym_fixed < 1e-6) && isTRUE(asym_buggy > 1e-3) &&
      isTRUE(ll_fixed <= ll_buggy)
cat("\n")
if (!ok) {
  cat("VERDICT: unexpected -- inspect (the bug should break antisymmetry while the fix preserves it).\n")
  quit(status = 1L)
}
cat("VERDICT: the buggy ancestor violates label-flip antisymmetry (provably wrong);\n",
    "        the REG-01 fix restores it and lowers held-out logistic loss.\n", sep = "")
