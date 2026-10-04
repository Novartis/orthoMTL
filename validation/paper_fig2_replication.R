#!/usr/bin/env Rscript
# =====================================================================
# Replication of Fig. 2 of Vervier et al. (2014), "On Learning Matrices
# with Orthogonal Columns or Disjoint Supports" (ECML-PKDD), with the
# CURRENT orthoMTL solver.
# ---------------------------------------------------------------------
# Source of the design: orthopen/vignettes/orthogonality_example.md
# (section "Convexity impact on model columns orthogonality") and paper
# Sec. 6.1. Question: how does the convexity of the penalty -- controlled
# by the diagonal of K = (off-diagonal 1, diagonal `delta`) -- affect
# (a) test MSE and (b) how orthogonal the learned columns are?
#
# Convexity (paper eq. 20): delta > T-1 strictly convex; delta = T-1 the
# atomic-norm limit (dotted line in the figure); delta < T-1 NON-convex.
#
# Setting: T = 10 tasks, d = 10 features, n_train = 50, n_test = 1000,
# W = random matrix with orthogonal columns, Y_train = X W + noise*N(0,1),
# Y_test noiseless. For each (noise, delta) the best lambda is picked by
# TEST MSE (as in the paper), repeated `reps` times.
#
# Deviations from the original vignette (which does not run as written):
#   * results indexed by noise level (the original indexed a list with
#     2.5 and re-used the loop variable after the loop);
#   * orthoMTL() instead of orthopen() (`$B`; `$W` is a back-compat alias);
#   * NaN correlations (all-zero columns at very large lambda) are treated
#     as "no deviation" instead of poisoning the average (such lambdas are
#     never the test-MSE argmin anyway);
#   * `reps` and the schedule are parameters (paper: 100 repeats).
# NOISE CONVENTION: the paper text says noise "variance" in {1, 2.5, 4}
# but the vignette (which generated the figure) uses `noise` as a
# multiplier of N(0,1), i.e. a standard deviation. The paper's reported
# test-MSE ranges (~0.6-0.9, ~2.75-3.0, ~3.8-4.1) are only consistent with
# the sd convention, so that is the default here.
#
# Run (against the INSTALLED package):
#   Rscript validation/paper_fig2_replication.R
# Env: FIG2_REPS (default 20), FIG2_CORES (default 4),
#      FIG2_SCHEDULE (default "sqrt": faithful; non-convex fits can land on
#      different stationary points under other schedules), FIG2_OUT (dir
#      for the PNG + RDS, default tempdir()).
# =====================================================================

suppressMessages(library(orthoMTL))

reps     <- as.integer(Sys.getenv("FIG2_REPS", "20"))
n_cores  <- as.integer(Sys.getenv("FIG2_CORES", "4"))
sched    <- Sys.getenv("FIG2_SCHEDULE", "sqrt")
out_dir  <- Sys.getenv("FIG2_OUT", tempdir())

Tt      <- 10L; NTRAIN <- 50L; NVAR <- 10L; NTEST <- 1000L
step    <- 0.1
DIAGVAL <- seq(1, 2 * Tt, by = 2)
LAMBDAS <- 1.6^seq(-16, 4)
NOISES  <- c(1, 2.5, 4)

one_rep <- function(task) {
  noise <- task$noise; irep <- task$irep
  suppressMessages(library(orthoMTL))
  set.seed(100000 * match(noise, NOISES) + irep)
  Xtr <- matrix(rnorm(NTRAIN * NVAR), NTRAIN, NVAR)
  Xte <- matrix(rnorm(NTEST  * NVAR), NTEST,  NVAR)
  W   <- qr.Q(qr(matrix(rnorm(NVAR * Tt), NVAR, Tt)))
  Ytr <- Xtr %*% W + matrix(rnorm(NTRAIN * Tt), NTRAIN) * noise
  Yte <- Xte %*% W
  K   <- matrix(1, Tt, Tt)

  mse <- matrix(NA_real_, length(LAMBDAS), length(DIAGVAL))
  ang <- matrix(NA_real_, length(LAMBDAS), length(DIAGVAL))
  for (id in seq_along(DIAGVAL)) {
    diag(K) <- DIAGVAL[id]
    for (il in seq_along(LAMBDAS)) {
      fit <- tryCatch(
        orthoMTL(Xtr, Ytr, lambda = LAMBDAS[il], K = K, step_size = step,
                 max_iter = 1e5, schedule = sched, verbose = 0),
        error = function(e) NULL)
      if (is.null(fit)) next
      B <- fit$B
      mse[il, id] <- 0.5 * sum((Xte %*% B - Yte)^2) / NTEST
      v <- suppressWarnings(acos(pmin(pmax(cor(B), -1), 1))) * 180 / pi - 90
      v[is.na(v)] <- 0
      diag(v) <- 0
      ang[il, id] <- sum(abs(v)) / (Tt * (Tt - 1))
    }
  }
  best <- apply(mse, 2, which.min)               # best lambda per diagonal
  list(noise = noise, irep = irep,
       mse = vapply(seq_along(DIAGVAL), function(i) mse[best[i], i], numeric(1)),
       ang = vapply(seq_along(DIAGVAL), function(i) ang[best[i], i], numeric(1)),
       lam = LAMBDAS[best])
}

tasks <- expand.grid(noise = NOISES, irep = seq_len(reps))
tasks <- lapply(seq_len(nrow(tasks)), function(i) as.list(tasks[i, ]))

t0 <- Sys.time()
rds_path <- file.path(out_dir, "fig2_replication.rds")
if (Sys.getenv("FIG2_REUSE", "0") == "1" && file.exists(rds_path)) {
  res <- readRDS(rds_path)$res            # re-plot/re-summarise without refitting
} else if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(p) .libPaths(p), .libPaths())
  parallel::clusterExport(cl, c("reps", "Tt", "NTRAIN", "NVAR", "NTEST", "step",
                                "DIAGVAL", "LAMBDAS", "NOISES", "sched"))
  res <- parallel::parLapplyLB(cl, tasks, one_rep)
} else {
  res <- lapply(tasks, one_rep)
}
cat(sprintf("done in %.1f min (reps=%d, schedule=%s)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), reps, sched))

# ---- aggregate --------------------------------------------------------
summ <- function(noise, field) {
  m <- do.call(rbind, lapply(Filter(function(r) r$noise == noise, res),
                             function(r) r[[field]]))
  ci <- 1.96 * apply(m, 2, sd, na.rm = TRUE) / sqrt(nrow(m))
  ci[is.na(ci)] <- 0                              # single repeat: no spread
  list(mean = colMeans(m, na.rm = TRUE), ci = ci)
}
# Paired comparison: every diagonal within a repeat sees the SAME data, so the
# difference to a reference diagonal removes between-dataset variance that
# inflates the unpaired error bars. Reference = the atomic-norm limit T-1.
REF <- which(DIAGVAL == Tt - 1)
paired <- function(noise) {
  m <- do.call(rbind, lapply(Filter(function(r) r$noise == noise, res),
                             function(r) r$mse))
  d <- m - m[, REF]
  ci <- 1.96 * apply(d, 2, sd, na.rm = TRUE) / sqrt(nrow(d))
  ci[is.na(ci)] <- 0
  list(mean = colMeans(d, na.rm = TRUE), ci = ci)
}
out <- list()
for (nz in NOISES) out[[as.character(nz)]] <- list(mse = summ(nz, "mse"),
                                                   ang = summ(nz, "ang"),
                                                   dmse = paired(nz))
saveRDS(list(res = res, summary = out, diag = DIAGVAL),
        file.path(out_dir, "fig2_replication.rds"))

# ---- figure (same layout as the paper: MSE top, angle bottom) ----------
png(file.path(out_dir, "fig2_replication.png"), width = 1500, height = 1150, res = 120)
op <- par(mfrow = c(3, 3), mar = c(4, 4.5, 2.5, 0.8))
ylabs <- c(mse = "Test MSE (unpaired 95% CI)",
           dmse = sprintf("Test MSE - MSE(diag=%d)  (paired 95%% CI)", Tt - 1),
           ang = "mean(|angle-90|)")
for (field in c("mse", "dmse", "ang")) for (nz in NOISES) {
  s <- out[[as.character(nz)]][[field]]
  yl <- range(c(s$mean + s$ci, s$mean - s$ci))
  plot(DIAGVAL, s$mean, type = "l", lwd = 3, ylim = yl, xlab = "Diagonal",
       ylab = ylabs[[field]], cex.lab = 0.9,
       main = if (field == "mse") sprintf("noise = %g", nz) else "")
  arrows(DIAGVAL, s$mean - s$ci, DIAGVAL, s$mean + s$ci, angle = 90, code = 3,
         length = 0.03)
  abline(v = Tt - 1, lty = 2, lwd = 2); grid()
  if (field == "dmse") abline(h = 0, col = "grey40")
}
par(op); dev.off()

# ---- report ------------------------------------------------------------
cat(sprintf("\nFig. 2 replication  (T=%d d=%d n=%d, %d repeats, noise as sd)\n",
            Tt, NVAR, NTRAIN, reps))
cat("Dotted line in the paper at diagonal = T-1 = 9 (left: non-convex).\n\n")
for (nz in NOISES) {
  s <- out[[as.character(nz)]]
  cat(sprintf("noise %g:\n  diag  :", nz), sprintf("%6g", DIAGVAL), "\n")
  cat("  MSE   :", sprintf("%6.3f", s$mse$mean), "\n")
  cat("  angle :", sprintf("%6.2f", s$ang$mean), "\n")
  cat(sprintf("  best diagonal (min MSE) = %g\n", DIAGVAL[which.min(s$mse$mean)]))
  cat("  paired dMSE vs diag", Tt - 1, ":", sprintf("%+.3f", s$dmse$mean), "\n")
  cat("  paired 95% CI half-w :", sprintf("%6.3f", s$dmse$ci), "\n")
  sig <- ifelse(s$dmse$mean + s$dmse$ci < 0, "better",
                ifelse(s$dmse$mean - s$dmse$ci > 0, "worse", "n.s."))
  cat("  vs diag", Tt - 1, "        :", sprintf("%6s", sig), "\n\n")
}

# ---- paper's qualitative claims, checked -------------------------------
best_d <- vapply(NOISES, function(nz)
  DIAGVAL[which.min(out[[as.character(nz)]]$mse$mean)], numeric(1))
ang_ok <- vapply(NOISES, function(nz) {
  a <- out[[as.character(nz)]]$ang$mean
  mean(a[DIAGVAL < Tt - 1]) < mean(a[DIAGVAL > Tt - 1])
}, logical(1))
d1 <- out[["1"]]$dmse
claims <- c(
  "low noise: best diagonal is non-convex (< T-1)"            = best_d[1] < Tt - 1,
  "low noise: diag 1 beats the atomic limit, PAIRED CI excludes 0" =
    d1$mean[1] + d1$ci[1] < 0,
  "best diagonal moves up (more convex) as noise grows"        = best_d[3] >= best_d[1],
  "non-convex K gives more orthogonal columns (all 3 noises)"  = all(ang_ok)
)
cat("Paper claims reproduced?\n")
for (nm in names(claims)) cat(sprintf("  [%s] %s\n", if (claims[[nm]]) "yes" else "NO ", nm))
cat(sprintf("\nFigure + data: %s\n", out_dir))
if (!all(claims)) quit(status = 1L)
