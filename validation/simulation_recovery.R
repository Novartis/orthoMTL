#!/usr/bin/env Rscript
# =====================================================================
# Tier-B confidence: survival-path ground-truth recovery study
# ---------------------------------------------------------------------
# The parity check (validation/parity_orthopen.R) proves orthoMTL's core
# solver is bit-identical to the validated orthopen v1.1.0 -- but only on
# the cases the two packages share. orthoMTL's SURVIVAL masking
# (survival = TRUE, censored.mat) is orthoMTL-only and parity cannot reach
# it. This study fills that gap.
#
# NOTE ON SCOPE: the MSE survival path exercised here is UNCHANGED by the
# REG-01/02/03 fixes (those touched the logistic loss/gradient and the
# alpha penalty; alpha = 0 reproduces the pre-fix default exactly). So this
# study is not an "old vs fixed" A/B -- the logistic A/B is covered by the
# orthopen parity check. Its job is to confirm the survival path EXTRACTS
# REAL, CORRECTLY-ORIENTED SIGNAL and behaves sanely, using the known
# ground truth from simulate_mtl().
#
# Metrics, each paired real-vs-null on the same data/split (null = outcome
# permuted across patients, breaking the X<->outcome link):
#   * support recovery -- AUC of |B_hat| row-norm vs the true signal set
#                         (the primary, most robust signal)
#   * held-out C-index -- cindex_mtl() on an untouched test split
#
# CEILING CONTEXT: we also report the C-index of the TRUE generative
# coefficients ("hazard-oracle"). It is LOW/anti-concordant by design --
# simulate_mtl() puts effects on the hazard scale (positive beta -> more
# hazard -> SHORTER survival) while the labels encode progression-free
# (1 = survived past threshold), and cindex_mtl() row-sums across tasks,
# collapsing the time-varying effect structure. The fitted model beating
# this oracle at label prediction is itself a correctness signal.
#
# Run:
#   "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/simulation_recovery.R
# Requires orthoMTL INSTALLED.
# =====================================================================

suppressMessages(library(orthoMTL))

# ---- study configuration (fixed, not tuned-to-win) ------------------
n          <- 200L
p          <- 30L
n_signals  <- 5L
thresholds <- c(4, 6, 10, 15)
lambda     <- 1e-3                 # the value used in the package example
train_frac <- 0.70
seeds      <- 1:10

# ---- helpers --------------------------------------------------------
# Mann-Whitney AUC (no extra package dependency)
auc <- function(score, pos) {
  pos <- as.logical(pos)
  np <- sum(pos); nn <- sum(!pos)
  if (np == 0L || nn == 0L) return(NA_real_)
  r <- rank(score)
  (sum(r[pos]) - np * (np + 1) / 2) / (np * nn)
}

# C-index of a fixed coefficient matrix B on (X, SurvTime, Event)
cindex_of_B <- function(B, X, S, E, thresholds) {
  Mt <- t(apply(X %*% B, 1, nnmaxheap_C))      # survival monotonicity projection
  cindex_mtl(create_longitudinal_labels(S, E, thresholds), Mt)
}

# Fit the documented survival workflow on a train split, score on test.
fit_and_score <- function(X, SurvTime, Event, thresholds, lambda,
                          train_idx, true_signal) {
  Xtr <- X[train_idx, , drop = FALSE]
  Ytr <- create_longitudinal_labels(SurvTime[train_idx], Event[train_idx], thresholds)
  Wtr <- create_indicator_matrix(Ytr)          # 1 = observed, 0 = censored
  K   <- create_constraint_matrix(length(thresholds))

  fit <- orthoMTL(Xtr, Ytr, lambda = lambda, K = K,
                  survival = TRUE, censored.mat = Wtr,
                  max_iter = 50000, verbose = 0)

  Yte   <- create_longitudinal_labels(SurvTime[-train_idx], Event[-train_idx], thresholds)
  preds <- predict(fit, newdata = X[-train_idx, , drop = FALSE])

  list(auc       = auc(rowSums(abs(fit$B)), true_signal),
       cindex    = cindex_mtl(Yte, preds),
       converged = isTRUE(fit$converged))
}

# ---- run over seeds: real model + permuted-outcome null -------------
rows <- vector("list", length(seeds))
for (k in seq_along(seeds)) {
  s   <- seeds[k]
  sim <- simulate_mtl(n = n, p = p, n_signals = n_signals,
                      thresholds = thresholds, seed = s)

  # true signal = any feature with a non-zero ground-truth coefficient
  true_signal <- rowSums(abs(sim$ground_truth$coefficients)) > 0

  set.seed(1000 + s)
  train_idx <- sort(sample(seq_len(n), floor(train_frac * n)))

  real <- fit_and_score(sim$X, sim$SurvTime, sim$Event,
                        thresholds, lambda, train_idx, true_signal)

  perm <- sample(seq_len(n))       # break X <-> outcome link
  null <- fit_and_score(sim$X, sim$SurvTime[perm], sim$Event[perm],
                        thresholds, lambda, train_idx, true_signal)

  # ceiling context: C-index of the true generative coefficients (test split)
  oracle <- cindex_of_B(sim$ground_truth$coefficients,
                        sim$X[-train_idx, , drop = FALSE],
                        sim$SurvTime[-train_idx], sim$Event[-train_idx], thresholds)

  rows[[k]] <- data.frame(
    seed = s,
    auc_real = real$auc,     auc_null = null$auc,
    cidx_real = real$cindex, cidx_null = null$cindex, cidx_oracle = oracle,
    converged = real$converged
  )
}
res <- do.call(rbind, rows)

# ---- report ---------------------------------------------------------
fmt <- function(x) sprintf("%.3f ± %.3f", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
n_pos <- sum(rowSums(abs(simulate_mtl(n = n, p = p, n_signals = n_signals,
                                      thresholds = thresholds, seed = seeds[1]
                                      )$ground_truth$coefficients)) > 0)

cat("Survival-path ground-truth recovery (orthoMTL, survival=TRUE, MSE loss)\n")
cat(sprintf("n=%d  p=%d  true-signal features=%d  thresholds=%s  lambda=%g\n",
            n, p, n_pos, paste(thresholds, collapse = ","), lambda))
cat(sprintf("%d seeds, %.0f/%.0f train/test split  (paired real-vs-null)\n\n",
            length(seeds), train_frac * 100, (1 - train_frac) * 100))

cat(sprintf("%-22s %-16s %-16s\n", "", "real model", "null (permuted)"))
cat(sprintf("%-22s %-16s %-16s   real>null in %d/%d seeds\n",
            "support-recovery AUC", fmt(res$auc_real), fmt(res$auc_null),
            sum(res$auc_real > res$auc_null), nrow(res)))
cat(sprintf("%-22s %-16s %-16s   real>null in %d/%d seeds\n",
            "held-out C-index", fmt(res$cidx_real), fmt(res$cidx_null),
            sum(res$cidx_real > res$cidx_null), nrow(res)))
cat(sprintf("%-22s %-16s   (low/anti-concordant by design -- see header)\n",
            "  hazard-oracle C-idx", fmt(res$cidx_oracle)))
cat(sprintf("%-22s %d/%d converged\n", "convergence", sum(res$converged), nrow(res)))

# ---- verdict --------------------------------------------------------
auc_wins  <- sum(res$auc_real  > res$auc_null)
cidx_wins <- sum(res$cidx_real > res$cidx_null)
beats_oracle <- mean(res$cidx_real, na.rm = TRUE) > mean(res$cidx_oracle, na.rm = TRUE)
cat(sprintf("\nsupport-AUC real>null: %d/%d   C-index real>null: %d/%d   model beats hazard-oracle: %s\n",
            auc_wins, nrow(res), cidx_wins, nrow(res), beats_oracle))

# Primary criterion: support recovery is the robust signal. C-index is a
# weaker secondary (low ceiling); we require it to lean the right way and
# the fitted model to out-predict the raw generative coefficients.
pass <- auc_wins >= 8L &&
        mean(res$auc_real, na.rm = TRUE) - mean(res$auc_null, na.rm = TRUE) > 0.10 &&
        cidx_wins >= 6L &&
        beats_oracle
if (!pass) {
  cat("VERDICT: signal weaker than expected -- inspect before trusting.\n")
  quit(status = 1L)
}
cat("VERDICT: survival path recovers real, correctly-oriented signal;",
    "support recovery is well above the null and the model out-predicts",
    "the generative coefficients on the label task.\n")
