#!/usr/bin/env Rscript
# =====================================================================
# SUPERSEDED by ab_k01b_K_shape_convexity_matched.R -- kept as a record.
# This first attempt is CONFOUNDED BY CONVEXITY: with T = 8 tasks the
# convexity threshold on diag(K) is T-1 = 7 (Vervier et al. 2014), yet the
# diagonal grid was {0.25, 1, 4}, so every non-identity K was non-convex
# while the identity baseline is convex -- at high noise, where convex
# penalties are expected to win. Its "identity wins" outcome is therefore
# NOT evidence that K shape is irrelevant.
# =====================================================================
# A/B (K-01): does the SHAPE of the task-relationship matrix K matter?
# ---------------------------------------------------------------------
# KANBAN K-01 (a feature requested by colleagues/discussants, not
# maintainer-originated): before building a `make_K()` helper / K-grid in
# cv_orthoMTL, check whether K's shape changes results at all. Today the
# default is K = I (pure per-task ridge, NO orthogonality pressure) and
# cv_orthoMTL only tunes K's DIAGONAL.
#
# Controlled regression data with KNOWN support structure (the shipped
# simulate_mtl() templates give every signal a dense effect across tasks,
# so it cannot distinguish K shapes):
#   sliding : each task's support is a window that slides with the task
#             index -> adjacent tasks overlap, distant tasks are disjoint
#             (the prior the diffusion K encodes)
#   blocks  : two groups of tasks; support shared within a group,
#             disjoint across groups
#   shared  : all tasks share ONE support (negative control: an
#             orthogonality prior is WRONG here, identity should do fine)
#
# K candidates (every one tuned fairly: lambda x diag_val grid via
# cv_orthoMTL, 3-fold, then refit on the full training set):
#   identity  : shipped default
#   uniform   : off-diagonal 1 everywhere (orthogonalise all pairs)
#   diffusion : create_constraint_matrix() (distance-weighted)
#   kernel    : 1 - exp(-|i-j| / 2) (smooth distance weighting)
#   oracle    : 1 where TRUE supports are disjoint, 0 where they overlap
#               (upper bound on what a perfect K could do; uses the truth)
#
# Metrics (held-out, n_test = 1000): RMSE and mean per-task support AUC
# (|B_hat[, t]| scoring the true non-zeros of task t).
#
# Run (against the INSTALLED package; ~minutes):
#   "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/ab_k01_K_shape.R
# Env: ABK_SEEDS (default 10) to shorten a smoke test.
# =====================================================================

suppressMessages(library(orthoMTL))

n_seeds   <- as.integer(Sys.getenv("ABK_SEEDS", "10"))
n_train   <- 100L
n_test    <- 1000L
p         <- 80L             # wide + weak signal so support recovery is NOT saturated
k         <- 8L
win       <- 5L              # support size per task
noise_sd  <- 1
lambdas   <- c(1e-3, 1e-2, 1e-1)
diag_vals <- c(0.25, 1, 4)
folds_k   <- 3L
sched     <- "log"           # same optimum, far fewer iterations (KANBAN S-02)

# ---- true coefficient matrices (p x k) -------------------------------
make_B <- function(scenario, rng_seed) {
  set.seed(rng_seed)
  B <- matrix(0, p, k)
  sgn <- function(m) sample(c(-1, 1), m, replace = TRUE) * runif(m, 0.4, 0.6)
  if (scenario == "sliding") {
    for (t in seq_len(k)) {
      idx <- seq(1 + 3 * (t - 1), length.out = win)
      B[idx, t] <- sgn(win)
    }
  } else if (scenario == "blocks") {
    grp <- rep(1:2, each = k / 2)
    sup <- list(1:win, (win + 1):(2 * win))
    for (t in seq_len(k)) B[sup[[grp[t]]], t] <- sgn(win)
  } else if (scenario == "shared") {
    sup <- 1:win
    for (t in seq_len(k)) B[sup, t] <- sgn(win)
  }
  B
}

# ---- K candidates -----------------------------------------------------
make_Ks <- function(Btrue) {
  d <- abs(outer(seq_len(k), seq_len(k), "-"))
  uni <- matrix(1, k, k)
  ker <- 1 - exp(-d / 2)
  S <- (Btrue != 0) * 1
  ov <- crossprod(S) > 0                      # TRUE where supports overlap
  orc <- matrix(as.numeric(!ov), k, k)
  list(identity  = diag(1, k),
       uniform   = uni,
       diffusion = create_constraint_matrix(k),
       kernel    = ker,
       oracle    = orc)
}

# ---- metrics ----------------------------------------------------------
auc <- function(score, pos) {
  pos <- as.logical(pos); np <- sum(pos); nn <- sum(!pos)
  if (np == 0L || nn == 0L) return(NA_real_)
  r <- rank(score)
  (sum(r[pos]) - np * (np + 1) / 2) / (np * nn)
}
support_auc <- function(Bhat, Btrue)
  mean(vapply(seq_len(k), function(t) auc(abs(Bhat[, t]), Btrue[, t] != 0),
              numeric(1)), na.rm = TRUE)

# ---- one (scenario, seed) run -----------------------------------------
run_one <- function(scenario, seed) {
  Btrue <- make_B(scenario, 7000 + seed)
  set.seed(seed)
  Xtr <- matrix(rnorm(n_train * p), n_train, p)
  Xte <- matrix(rnorm(n_test  * p), n_test,  p)
  Ytr <- Xtr %*% Btrue + matrix(rnorm(n_train * k, sd = noise_sd), n_train, k)
  Yte <- Xte %*% Btrue + matrix(rnorm(n_test  * k, sd = noise_sd), n_test,  k)
  folds <- sample(rep(seq_len(folds_k), length.out = n_train))

  Ks <- make_Ks(Btrue)
  out <- lapply(names(Ks), function(nm) {
    r <- tryCatch({
      cv <- cv_orthoMTL(Xtr, Ytr, K = Ks[[nm]], lambdas = lambdas,
                        alphas = 0, stepsizes = 0.1, diag_vals = diag_vals,
                        survival = FALSE, logistic = FALSE, metric = "rmse",
                        schedule = sched, folds = folds, n_cores = 1,
                        seed = seed, verbose = FALSE)
      Kb <- Ks[[nm]]; diag(Kb) <- cv$best$diag_val
      fit <- orthoMTL(Xtr, Ytr, lambda = cv$best$lambda, alpha = 0,
                      step_size = cv$best$stepsize, K = Kb,
                      survival = FALSE, logistic = FALSE, schedule = sched,
                      seed = seed, max_iter = 50000, verbose = 0)
      data.frame(rmse = sqrt(mean((Xte %*% fit$B - Yte)^2)),
                 sauc = support_auc(fit$B, Btrue),
                 lambda = cv$best$lambda, diag_val = cv$best$diag_val)
    }, error = function(e) {
      data.frame(rmse = NA_real_, sauc = NA_real_, lambda = NA_real_,
                 diag_val = NA_real_)
    })
    cbind(scenario = scenario, seed = seed, K = nm, r)
  })
  do.call(rbind, out)
}

# ---- run ---------------------------------------------------------------
scenarios <- c("sliding", "blocks", "shared")
t0 <- Sys.time()
res <- do.call(rbind, lapply(scenarios, function(sc)
  do.call(rbind, lapply(seq_len(n_seeds), function(s) {
    cat(sprintf("  scenario=%s seed=%d (%.0fs)\n", sc, s,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    run_one(sc, s)
  }))))

# ---- report ------------------------------------------------------------
cat("\nA/B K-01: does K shape matter?  (regression mode, known support)\n")
cat(sprintf("n_train=%d p=%d tasks=%d support/task=%d noise_sd=%g seeds=%d\n",
            n_train, p, k, win, noise_sd, n_seeds))
cat("Each K tuned over lambda x diag_val by 3-fold CV; held-out n=1000.\n\n")

fmt <- function(x) sprintf("%.3f ± %.3f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
for (sc in scenarios) {
  d <- res[res$scenario == sc, ]
  cat(sprintf("== %s ==\n", sc))
  cat(sprintf("%-10s %-16s %-16s %-14s %s\n", "K", "test RMSE", "support AUC",
              "paired dRMSE vs identity (mean ± SE, wins)", ""))
  base <- d[d$K == "identity", ]
  for (nm in unique(d$K)) {
    e <- d[d$K == nm, ]
    diffs <- e$rmse - base$rmse[match(e$seed, base$seed)]
    se <- sd(diffs, na.rm = TRUE) / sqrt(sum(!is.na(diffs)))
    cat(sprintf("%-10s %-16s %-16s %+.3f ± %.3f  (better in %d/%d)\n",
                nm, fmt(e$rmse), fmt(e$sauc), mean(diffs, na.rm = TRUE), se,
                sum(diffs < 0, na.rm = TRUE), sum(!is.na(diffs))))
  }
  cat("\n")
}
if (any(is.na(res$rmse))) cat("NOTE: ", sum(is.na(res$rmse)), " failed fits (NA).\n")

# ---- verdict: is a non-identity, REALISTIC (non-oracle) K worth it? ---
sig <- function(sc) {
  d <- res[res$scenario == sc & res$K != "oracle", ]
  base <- d[d$K == "identity", ]
  best <- NULL
  for (nm in setdiff(unique(d$K), "identity")) {
    e <- d[d$K == nm, ]
    diffs <- e$rmse - base$rmse[match(e$seed, base$seed)]
    m <- mean(diffs, na.rm = TRUE)
    se <- sd(diffs, na.rm = TRUE) / sqrt(sum(!is.na(diffs)))
    if (is.null(best) || m < best$m) best <- list(K = nm, m = m, se = se)
  }
  best
}
cat("Verdict (best realistic non-identity K vs identity, paired dRMSE; <0 = better):\n")
for (sc in scenarios) {
  b <- sig(sc)
  flag <- if (is.finite(b$se) && b$m < -2 * b$se) "K SHAPE MATTERS" else "no clear effect"
  cat(sprintf("  %-8s best=%-9s dRMSE=%+.3f ± %.3f  -> %s\n",
              sc, b$K, b$m, b$se, flag))
}
cat("\nInterpretation: if 'sliding'/'blocks' show a clear gain and 'shared'",
    "does not hurt, a K helper / K-grid is worth building (K-01).\n",
    "If all three say 'no clear effect', shrink K-01 to validate_K() + docs.\n")
saveRDS(res, file.path(tempdir(), "ab_k01_results.rds"))
