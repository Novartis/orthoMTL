#!/usr/bin/env Rscript
# =====================================================================
# A/B (K-01, redo): does the SHAPE of K matter once CONVEXITY is matched?
# ---------------------------------------------------------------------
# validation/ab_k01_K_shape.R compared K shapes at diagonal values
# {0.25, 1, 4} with T = 8 tasks. By the convexity condition of Vervier et
# al. (2014) -- see validation/paper_fig2_replication.R -- almost all of
# those were NON-convex, against a convex identity K, at high noise (the
# regime where the paper shows convex penalties win). That comparison was
# confounded by convexity. This script removes the confound.
#
# Convexity of  Omega_K(W)^2 = sum_i K_ii ||w_i||^2 + sum_{i!=j} K_ij |w_i'w_j|
# (K entrywise >= 0): sufficient condition = companion matrix
# Kbar = diag(K_ii) - offdiag(K_ij) is positive semidefinite. Writing
# K = O + d*I with O the off-diagonal part, Kbar = d*I - O, so
#     convex  <=>  d >= lambda_max(O).
# We therefore parametrise every shape by
#     rho = d / lambda_max(O)      (rho < 1 non-convex, = 1 atomic limit,
#                                    > 1 strictly convex)
# which is exactly the paper's delta / (T-1) for the uniform shape. O is
# rescaled to lambda_max(O) = T-1 so shapes share a scale.
#
# Shapes: uniform (the paper's K), diffusion (create_constraint_matrix),
# kernel (1 - exp(-|i-j|/2)), oracle (1 where TRUE supports are disjoint;
# upper bound, uses the truth; skipped when all supports overlap), plus the
# identity-K baseline (pure ridge, convex).
# Scenarios (known structure):
#   orthogonal : the paper's setting -- dense orthogonal W, T = p = 10, n = 50
#   sliding    : support window slides with the task (adjacent tasks overlap)
#   blocks     : two task groups with disjoint supports
#   shared     : one common support (negative control: orthogonality is wrong)
# For each (shape, rho) the best lambda is chosen by TEST MSE as in the paper
# (isolates the effect of K; optimistic equally for all configs). Metric: MSE
# against the noiseless test signal. Reported as the PAIRED difference to the
# identity baseline on the same data (negative = better than identity).
#
# NOTE: the "best rho" columns select among 8 rho values on the same data, so
# the Delta at best rho is slightly optimistic; the full curves are the honest
# picture.
#
# Run (against the INSTALLED package):
#   Rscript validation/ab_k01b_K_shape_convexity_matched.R
# Env: ABK_SEEDS (default 20), ABK_CORES (default 6), ABK_SCHEDULE (default
#      "sqrt": faithful for non-convex fits), ABK_OUT (results dir, default
#      validation/_results/k01b -- persistent, NOT tempdir).
#
# CHECKPOINTED / RESUMABLE: every (scenario, noise, seed) task is saved as
# its own file in ABK_OUT/tasks/ as soon as it finishes (atomic write).
# Re-running the script skips finished tasks, so a crash loses at most the
# in-flight tasks. Tasks are ordered seed-major so a partial run covers every
# scenario/noise cell evenly. The report aggregates whatever is on disk and
# prints how many seeds each cell has.
#   ABK_MAX_TASKS  run at most this many NEW tasks, then aggregate (run in
#                  small batches; default: all pending)
#   ABK_AGG_ONLY=1 skip fitting, just aggregate + plot what exists
# =====================================================================

suppressMessages(library(orthoMTL))

n_seeds  <- as.integer(Sys.getenv("ABK_SEEDS", "20"))
n_cores  <- as.integer(Sys.getenv("ABK_CORES", "6"))
sched    <- Sys.getenv("ABK_SCHEDULE", "sqrt")
out_dir  <- Sys.getenv("ABK_OUT", file.path("validation", "_results", "k01b"))
task_dir <- file.path(out_dir, "tasks")
dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
max_new  <- as.integer(Sys.getenv("ABK_MAX_TASKS", "0"))   # 0 = all pending
agg_only <- Sys.getenv("ABK_AGG_ONLY", "0") == "1"
task_file <- function(scen, noise, seed)
  file.path(task_dir, sprintf("%s__noise%s__seed%03d.rds", scen, noise, seed))

NTEST   <- 1000L
step    <- 0.1
RHOS    <- c(0.05, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2.5)
LAMBDAS <- 1.6^seq(-16, 6, by = 2)

SCEN <- list(
  orthogonal = list(T = 10L, p = 10L, n = 50L,  noise = c(1, 2.5, 4)),
  sliding    = list(T = 8L,  p = 60L, n = 100L, noise = c(0.5, 1, 2)),
  blocks     = list(T = 8L,  p = 60L, n = 100L, noise = c(0.5, 1, 2)),
  shared     = list(T = 8L,  p = 60L, n = 100L, noise = c(0.5, 1, 2))
)

make_B <- function(scen, cfg) {
  Tt <- cfg$T; p <- cfg$p; win <- 5L
  mag <- function(m) sample(c(-1, 1), m, replace = TRUE) * runif(m, 0.8, 1.2)
  B <- matrix(0, p, Tt)
  if (scen == "orthogonal") {
    B <- qr.Q(qr(matrix(rnorm(p * Tt), p, Tt)))
  } else if (scen == "sliding") {
    for (t in seq_len(Tt)) B[seq(1 + 3 * (t - 1), length.out = win), t] <- mag(win)
  } else if (scen == "blocks") {
    grp <- rep(1:2, each = Tt / 2); sup <- list(1:win, (win + 1):(2 * win))
    for (t in seq_len(Tt)) B[sup[[grp[t]]], t] <- mag(win)
  } else if (scen == "shared") {
    for (t in seq_len(Tt)) B[1:win, t] <- mag(win)
  }
  B
}

# off-diagonal parts O (zero diagonal), rescaled to lambda_max = T-1
shapes_O <- function(Btrue) {
  Tt <- ncol(Btrue)
  d <- abs(outer(seq_len(Tt), seq_len(Tt), "-"))
  S <- (Btrue != 0) * 1
  disjoint <- crossprod(S) == 0
  Os <- list(
    uniform   = matrix(1, Tt, Tt),
    diffusion = create_constraint_matrix(Tt, diag_val = 0),
    kernel    = 1 - exp(-d / 2),
    oracle    = matrix(as.numeric(disjoint), Tt, Tt)
  )
  Os <- lapply(Os, function(O) { diag(O) <- 0; O })
  Os <- Filter(function(O) max(abs(O)) > 0, Os)       # drop empty oracle
  lapply(Os, function(O) {
    lmax <- max(eigen(O, symmetric = TRUE, only.values = TRUE)$values)
    O * ((ncol(O) - 1) / lmax)
  })
}

best_lambda_mse <- function(Xtr, Ytr, Xte, Yte, K) {
  m <- vapply(LAMBDAS, function(l) {
    fit <- tryCatch(orthoMTL(Xtr, Ytr, lambda = l, K = K, step_size = step,
                             max_iter = 1e5, schedule = sched, verbose = 0),
                    error = function(e) NULL)
    if (is.null(fit)) NA_real_ else mean((Xte %*% fit$B - Yte)^2)
  }, numeric(1))
  min(m, na.rm = TRUE)
}

one_task <- function(task) {
  suppressMessages(library(orthoMTL))
  scen <- task$scenario; noise <- task$noise; seed <- task$seed
  f <- task_file(scen, noise, seed)
  if (file.exists(f)) return(invisible(f))                # already done
  cfg <- SCEN[[scen]]
  set.seed(1e6 * match(scen, names(SCEN)) + 1e4 * match(noise, cfg$noise) + seed)
  Btrue <- make_B(scen, cfg)
  Xtr <- matrix(rnorm(cfg$n * cfg$p), cfg$n, cfg$p)
  Xte <- matrix(rnorm(NTEST * cfg$p), NTEST, cfg$p)
  Ytr <- Xtr %*% Btrue + matrix(rnorm(cfg$n * cfg$T), cfg$n) * noise
  Yte <- Xte %*% Btrue

  rows <- list()
  rows[[1]] <- data.frame(shape = "identity", rho = NA_real_,
                          mse = best_lambda_mse(Xtr, Ytr, Xte, Yte, diag(1, cfg$T)))
  Os <- shapes_O(Btrue)
  for (nm in names(Os)) for (rho in RHOS) {
    K <- Os[[nm]]; diag(K) <- rho * (cfg$T - 1)         # lambda_max(O) = T-1
    rows[[length(rows) + 1]] <- data.frame(
      shape = nm, rho = rho,
      mse = best_lambda_mse(Xtr, Ytr, Xte, Yte, K))
  }
  res_task <- cbind(scenario = scen, noise = noise, seed = seed,
                    do.call(rbind, rows))
  tmp <- paste0(f, ".tmp", Sys.getpid())                  # atomic write
  saveRDS(res_task, tmp)
  file.rename(tmp, f)
  invisible(f)
}

grid  <- do.call(rbind, lapply(names(SCEN), function(s)
  expand.grid(scenario = s, noise = SCEN[[s]]$noise, seed = seq_len(n_seeds),
              stringsAsFactors = FALSE)))
grid  <- grid[order(grid$seed), ]                  # seed-major: even partial coverage
tasks <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))
done  <- vapply(tasks, function(t) file.exists(task_file(t$scenario, t$noise, t$seed)),
                logical(1))
pending <- tasks[!done]
if (max_new > 0L && length(pending) > max_new) pending <- pending[seq_len(max_new)]
cat(sprintf("tasks: %d total, %d already on disk, %d to run now (dir: %s)\n",
            length(tasks), sum(done), if (agg_only) 0L else length(pending), out_dir))

t0 <- Sys.time()
if (!agg_only && length(pending) > 0L) {
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(p) .libPaths(p), .libPaths())
    parallel::clusterExport(cl, c("NTEST", "step", "RHOS", "LAMBDAS", "SCEN",
                                  "sched", "make_B", "shapes_O",
                                  "best_lambda_mse", "task_file", "task_dir"))
    parallel::parLapplyLB(cl, pending, one_task)
  } else {
    lapply(pending, one_task)
  }
  cat(sprintf("batch done in %.1f min (schedule=%s)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), sched))
}

# ---- load every finished task from disk (partial results are fine) -----
files <- list.files(task_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0L) stop("no finished tasks in ", task_dir, call. = FALSE)
res <- do.call(rbind, lapply(files, readRDS))
seeds_per_cell <- aggregate(seed ~ scenario + noise, res[res$shape == "identity", ],
                            function(s) length(unique(s)))
cat("seeds available per cell:\n"); print(seeds_per_cell, row.names = FALSE)
cat(sprintf("(target: %d seeds per cell; %d of %d tasks finished)\n",
            n_seeds, length(files), length(tasks)))

# ---- paired difference vs identity (same scenario/noise/seed) ----------
base <- res[res$shape == "identity", c("scenario", "noise", "seed", "mse")]
names(base)[4] <- "mse_id"
d <- merge(res[res$shape != "identity", ], base, by = c("scenario", "noise", "seed"))
d$dmse <- d$mse - d$mse_id

agg <- do.call(rbind, lapply(
  split(d, list(d$scenario, d$noise, d$shape, d$rho), drop = TRUE), function(g)
    data.frame(scenario = g$scenario[1], noise = g$noise[1], shape = g$shape[1],
               rho = g$rho[1], dmse = mean(g$dmse),
               ci = 1.96 * sd(g$dmse) / sqrt(nrow(g)), n = nrow(g),
               mse = mean(g$mse), mse_id = mean(g$mse_id))))
agg$ci[is.na(agg$ci)] <- 0                      # single seed: no spread
agg$scenario <- factor(agg$scenario, names(SCEN))
# noise columns by LEVEL (the orthogonal scenario uses a different sd set)
agg$level <- factor(mapply(function(s, nz) match(nz, SCEN[[as.character(s)]]$noise),
                           agg$scenario, agg$noise),
                    levels = 1:3, labels = c("low noise", "mid noise", "high noise"))

# ---- figure -------------------------------------------------------------
suppressMessages(library(ggplot2))
p <- ggplot(agg, aes(rho, dmse, colour = shape, fill = shape)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_ribbon(aes(ymin = dmse - ci, ymax = dmse + ci), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  scale_x_log10(breaks = RHOS, labels = RHOS) +
  facet_grid(scenario ~ level, scales = "free_y") +
  labs(x = expression(rho == diag(K) / lambda[max](O) ~ " (<1 non-convex, 1 = atomic limit, >1 convex)"),
       y = "Test-MSE difference vs identity K (paired 95% CI; < 0 = better)",
       caption = paste0("Noise sd (low/mid/high): orthogonal 1 / 2.5 / 4; ",
                        "sliding, blocks, shared 0.5 / 1 / 2.  Best lambda per point ",
                        "chosen on test MSE (as in the paper).")) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom",
                                   axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "ab_k01b_K_shape.png"), p, width = 11, height = 10, dpi = 120)

# ---- report -------------------------------------------------------------
cat("\nPaired Delta test-MSE vs identity K (negative = better). For each shape:\n")
cat("best rho (lowest mean Delta), Delta there +- paired 95% CI, and whether that\n")
cat("best rho is non-convex (<1), atomic (=1) or convex (>1).\n\n")
for (sc in names(SCEN)) for (nz in SCEN[[sc]]$noise) {
  a <- agg[agg$scenario == sc & agg$noise == nz, ]
  cat(sprintf("%-10s noise %-4g  identity MSE = %.3f\n", sc, nz, a$mse_id[1]))
  for (sh in unique(a$shape)) {
    z <- a[a$shape == sh, ]; b <- z[which.min(z$dmse), ]
    verdict <- if (b$n < 3) sprintf("n=%d, too few seeds", b$n)
               else if (b$dmse + b$ci < 0) "BETTER"
               else if (b$dmse - b$ci > 0) "worse" else "n.s."
    cat(sprintf("   %-10s best rho=%-5g %-10s Delta=%+.4f +- %.4f  (%s)\n",
                sh, b$rho, if (b$rho < 1) "non-convex" else if (b$rho == 1) "atomic" else "convex",
                b$dmse, b$ci, verdict))
  }
}
cat(sprintf("\nFigure + data: %s\n", out_dir))
