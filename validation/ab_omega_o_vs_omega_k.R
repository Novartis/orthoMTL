#!/usr/bin/env Rscript
# =====================================================================
# A/B (O-01): the pairwise penalty Omega_K (Vervier et al. 2014, shipped
# default) versus the atomic norm Omega_O (penalty = "omega_O"), plus a lasso
# floor. What does each bring, where do they overlap, where do they differ?
# ---------------------------------------------------------------------
# Background (derivations in ai-research audit/ortho_penalty_screening,
# Proposition C and the task123 record). With G = W^T W:
#   Omega_K(W)^2 = sum_st K_st |G_st| = max over sign matrices S of tr(W (K o S) W^T)
#   Omega_O(W)^2 = max over correlation matrices Z (Z psd, diag 1) of tr(W Z W^T)
# Omega_K is pairwise and carries a prior K; it is convex only when the
# companion matrix of K is psd, which forces a ridge-like diagonal.
# Omega_O has no K, is convex for every T, and equals Omega_K with an all-ones
# K exactly when G is balanced; it is smaller on frustrated G.
#
# Parts (ABO_PARTS, comma-separated; default G,R,S,F; OG is opt-in and needs the gamma build):
#   G  penalty geometry: exact penalty values on canonical W (no learning).
#   R  K-01c regression scenarios (orthogonal / sliding / blocks / shared x 3
#      noise x seeds), same generators, seeds, folds and noiseless test set as
#      validation/ab_k01c_K_cv_realistic.R. Adds arms O (omega_O) and L (lasso),
#      lambda by the same 3-fold inner CV, and PAIRS them with the K-01c arms
#      stored on disk (A identity, B_* fixed shapes, C full CV over K,
#      D estimated K, E oracle K). Each task first refits K-01c's arm A and
#      asserts the test MSE is reproduced exactly (same data, same solver).
#   S  survival, vignette protocol (simulate_mtl, diffusion K, the vignette's
#      cv_orthoMTL grid with lambda extended to 1e-5..1e-1, 5 folds, C-index),
#      T = 4 and T = 8, held-out test half.
#   F  new regression scenarios probing the fundamentals: frustrated (two
#      triples at 120 degrees), rotating (temporal drift through a sign
#      switch), orthogonal20 (T = 20). Arms A, B_uniform, B_diffusion, C, O, L.
#
# Run (against a library holding the omega-o build; K-01c results needed for R):
#   R_LIBS=<lib> Rscript validation/ab_omega_o_vs_omega_k.R
# Env: ABO_PARTS, ABO_SEEDS (12), ABO_CORES (6), ABO_SCHEDULE ("sqrt", as K-01c),
#      ABO_OUT (validation/_results/omega_o), ABO_K01C (K-01c tasks dir,
#      validation/_results/k01c/tasks), ABO_AGG_ONLY=1.
# CHECKPOINTED / RESUMABLE: one atomic .rds per task under ABO_OUT/<part>/.
# =====================================================================

suppressMessages(library(orthoMTL))
stopifnot("needs the omega-o build (penalty= argument)" =
            "penalty" %in% names(formals(orthoMTL)))
cat(sprintf("orthoMTL %s from %s\n", as.character(packageVersion("orthoMTL")),
            find.package("orthoMTL")))

parts    <- strsplit(Sys.getenv("ABO_PARTS", "G,R,S,F"), ",")[[1]]
n_seeds  <- as.integer(Sys.getenv("ABO_SEEDS", "12"))
n_cores  <- as.integer(Sys.getenv("ABO_CORES", "6"))
sched    <- Sys.getenv("ABO_SCHEDULE", "sqrt")
out_dir  <- Sys.getenv("ABO_OUT", file.path("validation", "_results", "omega_o"))
k01c_dir <- Sys.getenv("ABO_K01C", file.path("validation", "_results", "k01c", "tasks"))
agg_only <- Sys.getenv("ABO_AGG_ONLY", "0") == "1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- shared constants and helpers (R and F copy K-01c exactly) ------------
NTEST   <- 1000L
step    <- 0.1
CV_K    <- 3L
RHOS    <- c(0.05, 0.4, 1, 1.5, 2.5)
LAMBDAS <- 1.6^seq(-14, 7, by = 3)                     # 8 values, as K-01c

SCEN_R <- list(
  orthogonal = list(T = 10L, p = 10L, n = 50L,  noise = c(1, 2.5, 4)),
  sliding    = list(T = 8L,  p = 60L, n = 100L, noise = c(0.5, 1, 2)),
  blocks     = list(T = 8L,  p = 60L, n = 100L, noise = c(0.5, 1, 2)),
  shared     = list(T = 8L,  p = 60L, n = 100L, noise = c(0.5, 1, 2))
)
SCEN_F <- list(
  frustrated   = list(T = 6L,  p = 30L, n = 100L, noise = c(0.5, 1, 2)),
  rotating     = list(T = 8L,  p = 40L, n = 100L, noise = c(0.5, 1, 2)),
  orthogonal20 = list(T = 20L, p = 20L, n = 100L, noise = c(1, 2.5, 4))
)

make_B <- function(scen, cfg) {                         # K-01c generators, verbatim
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

# new generators (part F). Each column has unit-ish norm, like K-01c's.
make_B_F <- function(scen, cfg) {
  Tt <- cfg$T; p <- cfg$p; win <- 5L
  plane <- function(rows) {                             # orthonormal pair on `rows`
    Q <- qr.Q(qr(matrix(rnorm(length(rows) * 2), length(rows), 2)))
    b1 <- b2 <- numeric(p); b1[rows] <- Q[, 1]; b2[rows] <- Q[, 2]
    list(b1 = b1, b2 = b2)
  }
  B <- matrix(0, p, Tt)
  if (scen == "frustrated") {
    # two groups of three tasks; inside a group the three coefficient vectors
    # sit at 120 degrees in a 2-D plane: every pairwise inner product is
    # negative (-1/2), the signed graph is frustrated (no flip makes all >= 0)
    for (g in 1:2) {
      pl <- plane(((g - 1) * 2 * win + 1):(g * 2 * win))
      ang <- c(0, 2, 4) * pi / 3 + runif(1, 0, 2 * pi)
      for (k in 1:3) B[, (g - 1) * 3 + k] <- cos(ang[k]) * pl$b1 + sin(ang[k]) * pl$b2
    }
  } else if (scen == "rotating") {
    # temporal drift: task t points at angle (t-1) * 150 / (T-1) degrees in one
    # plane, so neighbours are similar and the ends have opposite-sign overlap
    pl <- plane(1:(2 * win))
    ang <- (seq_len(Tt) - 1) * (150 / (Tt - 1)) * pi / 180
    for (t in seq_len(Tt)) B[, t] <- cos(ang[t]) * pl$b1 + sin(ang[t]) * pl$b2
  } else if (scen == "orthogonal20") {
    B <- qr.Q(qr(matrix(rnorm(p * Tt), p, Tt)))
  }
  B
}

scale_O <- function(O) {                                # K-01c, verbatim
  O[!is.finite(O)] <- 0; diag(O) <- 0
  if (max(abs(O)) == 0) return(NULL)
  lmax <- max(eigen(O, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(lmax) || lmax < 1e-8) return(NULL)
  O * ((ncol(O) - 1) / lmax)
}
K_from <- function(O, rho, Tt) { K <- O; diag(K) <- rho * (Tt - 1); K }

# one fit, scored as K-01c does; arm = list(penalty, K, alpha[, gamma])
# gamma is passed only when the arm sets it, so parts R/F/S keep running on a pre-gamma build
fit_mse <- function(Xtr, Ytr, Xev, Yev, lambda, arm) {
  args <- list(Xtr, Ytr, lambda = lambda, K = arm$K, alpha = arm$alpha,
               penalty = arm$penalty, step_size = step,
               max_iter = 1e5, schedule = sched, verbose = 0)
  if (!is.null(arm$gamma)) args$gamma <- arm$gamma
  fit <- tryCatch(do.call(orthoMTL, args), error = function(e) NULL)
  if (is.null(fit)) NA_real_ else mean((Xev %*% fit$B - Yev)^2)
}

# 3-fold CV over LAMBDAS for one arm, refit on all training, score once on test
cv_arm <- function(Xtr, Ytr, Xte, Yte, folds, arm) {
  cv <- vapply(LAMBDAS, function(l) mean(vapply(seq_len(CV_K), function(k) {
    tr <- folds != k; va <- folds == k
    fit_mse(Xtr[tr, , drop = FALSE], Ytr[tr, , drop = FALSE],
            Xtr[va, , drop = FALSE], Ytr[va, , drop = FALSE], l, arm)
  }, numeric(1))), numeric(1))
  i <- which.min(cv)
  c(lambda = LAMBDAS[i], cv = cv[i], mse = fit_mse(Xtr, Ytr, Xte, Yte, LAMBDAS[i], arm))
}

save_atomic <- function(x, f) {
  tmp <- paste0(f, ".tmp", Sys.getpid()); saveRDS(x, tmp); file.rename(tmp, f)
}

run_tasks <- function(tasks, fun, exports) {
  if (agg_only || length(tasks) == 0L) return(invisible(NULL))
  t0 <- Sys.time()
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(p) .libPaths(p), .libPaths())
    parallel::clusterExport(cl, exports, envir = globalenv())
    parallel::parLapplyLB(cl, tasks, fun)
  } else lapply(tasks, fun)
  cat(sprintf("  %d tasks in %.1f min\n", length(tasks),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

paired <- function(res, base = "A", keys = c("scenario", "noise", "seed"), value = "mse") {
  b <- res[res$strategy == base, c(keys, value)]; names(b)[length(b)] <- "base"
  d <- merge(res[res$strategy != base, ], b, by = keys)
  d$delta <- d[[value]] - d$base
  do.call(rbind, lapply(split(d, list(d$scenario, d$noise, d$strategy), drop = TRUE),
    function(g) data.frame(scenario = g$scenario[1], noise = g$noise[1],
                           strategy = g$strategy[1], n = nrow(g), delta = mean(g$delta),
                           ci = if (nrow(g) > 1) 1.96 * sd(g$delta) / sqrt(nrow(g)) else 0,
                           harm = mean(g$delta > 0), worst = max(g$delta),
                           base = mean(g$base))))
}
verdict <- function(z, lower_is_better = TRUE) {
  s <- if (lower_is_better) 1 else -1
  if (s * (z$delta + s * z$ci) < 0) "BETTER" else if (s * (z$delta - s * z$ci) > 0) "worse" else "n.s."
}
exports_common <- c("NTEST", "step", "CV_K", "RHOS", "LAMBDAS", "SCEN_R", "SCEN_F",
                    "sched", "make_B", "make_B_F", "scale_O", "K_from", "fit_mse",
                    "cv_arm", "save_atomic", "out_dir", "k01c_dir")

# =====================================================================
# Part G -- penalty geometry (exact, no learning)
# =====================================================================
if ("G" %in% parts) {
  cat("\n=== Part G: penalty values on canonical coefficient matrices ===\n")
  sdp  <- function(W) orthoMTL:::.maxcut_sdp(crossprod(W), max_sweeps = 5000)$value
  mcut <- function(W) { Tt <- ncol(W); E <- as.matrix(expand.grid(rep(list(c(-1, 1)), Tt)))
    max(apply(E, 1, function(e) sum((W %*% e)^2))) }
  pairK <- function(W, K) sum(abs(crossprod(W)) * K)
  conv  <- function(K) { C <- -abs(K); diag(C) <- diag(K)       # companion matrix
    min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) >= -1e-10 }
  u <- function(...) { v <- c(...); v / sqrt(sum(v^2)) }
  ang3 <- c(0, 2, 4) * pi / 3
  cfgs <- list(
    "orthogonal (3 tasks)"        = diag(3),
    "shared, same sign"           = cbind(u(1, 0, 0), u(1, 0, 0), u(1, 0, 0)),
    "shared, flipped signs"       = cbind(u(1, 0, 0), u(-1, 0, 0), u(1, 0, 0)),
    "partial overlap (cos .5)"    = cbind(u(1, 0, 0), u(1, sqrt(3), 0), u(0, 1, 0)),
    "frustrated 120 deg"          = rbind(cos(ang3), sin(ang3), 0),
    "frustrated mild (cos -.2)"   = { G <- matrix(-0.2, 3, 3); diag(G) <- 1; chol(G) },
    "rotation 0-135 deg (4 tasks)"= sapply(c(0, 45, 90, 135) * pi / 180, function(a) c(cos(a), sin(a), 0))
  )
  tab <- do.call(rbind, lapply(names(cfgs), function(nm) {
    W <- cfgs[[nm]]; Tt <- ncol(W); fro <- sum(W^2)
    J <- matrix(1, Tt, Tt)
    Ounif <- J; diag(Ounif) <- 0; Kconv <- Ounif; diag(Kconv) <- Tt - 1   # rho = 1
    Dd <- create_constraint_matrix(Tt, diag_val = 0); lmD <- max(eigen(Dd, TRUE, only.values = TRUE)$values)
    Kd_v <- create_constraint_matrix(Tt)                                  # vignette diag 0.5
    Kd_c <- Dd; diag(Kd_c) <- lmD                                          # convex boundary
    balanced <- isTRUE(all.equal(mcut(W), pairK(W, J)))
    data.frame(config = nm, T = Tt, balanced = balanced,
               ridge = fro, maxcut = mcut(W), omega_O = sdp(W),
               K_allones = pairK(W, J), K_unif_convex = pairK(W, Kconv),
               K_diff_vignette = pairK(W, Kd_v), K_diff_convex = pairK(W, Kd_c),
               check.names = FALSE)
  }))
  num <- vapply(tab, is.numeric, logical(1)) & names(tab) != "T"
  tab[num] <- lapply(tab[num], round, 3)
  print(tab, row.names = FALSE)
  Tt <- 3; J <- matrix(1, Tt, Tt); Kc <- J; diag(Kc) <- Tt - 1
  cat(sprintf(paste0("\nConvexity (companion-matrix test, T = 3): all-ones K %s; uniform at rho = 1 %s; ",
                     "vignette diffusion (diag 0.5) %s. omega_O: convex for every T (max of convex quadratics).\n"),
              if (conv(J)) "convex" else "NON-convex", if (conv(Kc)) "convex" else "NON-convex",
              if (conv(create_constraint_matrix(4))) "convex (T = 4)" else "NON-convex (T = 4)"))
  cat("Reading: ridge <= maxcut <= omega_O <= K_allones always; omega_O = K_allones iff balanced.\n")
  saveRDS(tab, file.path(out_dir, "part_G_geometry.rds"))
}

# =====================================================================
# Part R -- K-01c scenarios, arms O and L paired with the stored K arms
# =====================================================================
if ("R" %in% parts) {
  cat("\n=== Part R: K-01c regression scenarios (paired with stored K arms) ===\n")
  stopifnot("K-01c task results not found: set ABO_K01C" = dir.exists(k01c_dir))
  dR <- file.path(out_dir, "R"); dir.create(dR, showWarnings = FALSE)
  one_R <- function(task) {
    suppressMessages(library(orthoMTL))
    scen <- task$scenario; noise <- task$noise; seed <- task$seed
    f <- file.path(out_dir, "R", sprintf("%s__noise%s__seed%03d.rds", scen, noise, seed))
    if (file.exists(f)) return(invisible(f))
    kf <- file.path(k01c_dir, sprintf("%s__noise%s__seed%03d.rds", scen, noise, seed))
    stored <- readRDS(kf)
    cfg <- SCEN_R[[scen]]; Tt <- cfg$T
    # identical random stream to ab_k01c_K_cv_realistic.R
    set.seed(1e6 * match(scen, names(SCEN_R)) + 1e4 * match(noise, cfg$noise) + seed)
    Btrue <- make_B(scen, cfg)
    Xtr <- matrix(rnorm(cfg$n * cfg$p), cfg$n, cfg$p)
    Xte <- matrix(rnorm(NTEST * cfg$p), NTEST, cfg$p)
    Ytr <- Xtr %*% Btrue + matrix(rnorm(cfg$n * cfg$T), cfg$n) * noise
    Yte <- Xte %*% Btrue
    folds <- sample(rep(seq_len(CV_K), length.out = cfg$n))
    # pairing check: K-01c arm A (identity, its CV-chosen lambda) must reproduce
    a <- stored[stored$strategy == "A", ]
    mseA <- fit_mse(Xtr, Ytr, Xte, Yte, a$lambda,
                    list(penalty = "K", K = diag(1, Tt), alpha = 0))
    if (!isTRUE(all.equal(mseA, a$mse, tolerance = 1e-10)))
      stop(sprintf("pairing check failed for %s: refit A mse %.12g vs stored %.12g",
                   basename(f), mseA, a$mse))
    O <- cv_arm(Xtr, Ytr, Xte, Yte, folds, list(penalty = "omega_O", K = NULL, alpha = 0))
    L <- cv_arm(Xtr, Ytr, Xte, Yte, folds, list(penalty = "K", K = NULL, alpha = 1))
    out <- rbind(stored[, c("scenario", "noise", "seed", "strategy", "lambda", "cv", "mse")],
                 data.frame(scenario = scen, noise = noise, seed = seed,
                            strategy = c("O", "L"), lambda = c(O["lambda"], L["lambda"]),
                            cv = c(O["cv"], L["cv"]), mse = c(O["mse"], L["mse"])))
    save_atomic(out, f)
    invisible(f)
  }
  grid <- do.call(rbind, lapply(names(SCEN_R), function(s)
    expand.grid(scenario = s, noise = SCEN_R[[s]]$noise, seed = seq_len(n_seeds),
                stringsAsFactors = FALSE)))
  grid <- grid[order(grid$seed), ]
  tasks <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))
  pend <- Filter(function(t) !file.exists(file.path(dR, sprintf("%s__noise%s__seed%03d.rds",
                   t$scenario, t$noise, t$seed))), tasks)
  cat(sprintf("  %d tasks, %d pending\n", length(tasks), length(pend)))
  run_tasks(pend, one_R, c(exports_common, "one_R"))
}

# =====================================================================
# Part F -- new scenarios, all arms run here
# =====================================================================
if ("F" %in% parts) {
  cat("\n=== Part F: frustrated / rotating / orthogonal20 ===\n")
  dF <- file.path(out_dir, "F"); dir.create(dF, showWarnings = FALSE)
  one_F <- function(task) {
    suppressMessages(library(orthoMTL))
    scen <- task$scenario; noise <- task$noise; seed <- task$seed
    f <- file.path(out_dir, "F", sprintf("%s__noise%s__seed%03d.rds", scen, noise, seed))
    if (file.exists(f)) return(invisible(f))
    cfg <- SCEN_F[[scen]]; Tt <- cfg$T
    set.seed(2e6 + 1e5 * match(scen, names(SCEN_F)) + 1e4 * match(noise, cfg$noise) + seed)
    Btrue <- make_B_F(scen, cfg)
    Xtr <- matrix(rnorm(cfg$n * cfg$p), cfg$n, cfg$p)
    Xte <- matrix(rnorm(NTEST * cfg$p), NTEST, cfg$p)
    Ytr <- Xtr %*% Btrue + matrix(rnorm(cfg$n * cfg$T), cfg$n) * noise
    Yte <- Xte %*% Btrue
    folds <- sample(rep(seq_len(CV_K), length.out = cfg$n))
    shapes <- list(uniform = scale_O(matrix(1, Tt, Tt)),
                   diffusion = scale_O(create_constraint_matrix(Tt, diag_val = 0)))
    # K arms: (shape, rho, lambda) by the same inner CV as K-01c
    cand <- rbind(data.frame(shape = "identity", rho = NA_real_, lambda = LAMBDAS),
                  expand.grid(shape = names(shapes), rho = RHOS, lambda = LAMBDAS,
                              stringsAsFactors = FALSE))
    Kof <- function(s, r) if (s == "identity") diag(1, Tt) else K_from(shapes[[s]], r, Tt)
    cand$cv <- vapply(seq_len(nrow(cand)), function(i) mean(vapply(seq_len(CV_K), function(k) {
      tr <- folds != k; va <- folds == k
      fit_mse(Xtr[tr, , drop = FALSE], Ytr[tr, , drop = FALSE], Xtr[va, , drop = FALSE],
              Ytr[va, , drop = FALSE], cand$lambda[i],
              list(penalty = "K", K = Kof(cand$shape[i], cand$rho[i]), alpha = 0))
    }, numeric(1))), numeric(1))
    pick <- function(sel) { i <- which(sel & is.finite(cand$cv)); i[which.min(cand$cv[i])] }
    rowK <- function(st, i) data.frame(strategy = st, shape = cand$shape[i], rho = cand$rho[i],
      lambda = cand$lambda[i], cv = cand$cv[i],
      mse = fit_mse(Xtr, Ytr, Xte, Yte, cand$lambda[i],
                    list(penalty = "K", K = Kof(cand$shape[i], cand$rho[i]), alpha = 0)))
    rows <- rbind(rowK("A", pick(cand$shape == "identity")),
                  rowK("B_uniform", pick(cand$shape == "uniform")),
                  rowK("B_diffusion", pick(cand$shape == "diffusion")),
                  rowK("C", pick(rep(TRUE, nrow(cand)))))
    O <- cv_arm(Xtr, Ytr, Xte, Yte, folds, list(penalty = "omega_O", K = NULL, alpha = 0))
    L <- cv_arm(Xtr, Ytr, Xte, Yte, folds, list(penalty = "K", K = NULL, alpha = 1))
    rows <- rbind(rows, data.frame(strategy = c("O", "L"), shape = NA, rho = NA,
                                   lambda = c(O["lambda"], L["lambda"]), cv = c(O["cv"], L["cv"]),
                                   mse = c(O["mse"], L["mse"])))
    save_atomic(cbind(scenario = scen, noise = noise, seed = seed, rows), f)
    invisible(f)
  }
  grid <- do.call(rbind, lapply(names(SCEN_F), function(s)
    expand.grid(scenario = s, noise = SCEN_F[[s]]$noise, seed = seq_len(n_seeds),
                stringsAsFactors = FALSE)))
  grid <- grid[order(grid$seed), ]
  tasks <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))
  pend <- Filter(function(t) !file.exists(file.path(dF, sprintf("%s__noise%s__seed%03d.rds",
                   t$scenario, t$noise, t$seed))), tasks)
  cat(sprintf("  %d tasks, %d pending\n", length(tasks), length(pend)))
  run_tasks(pend, one_F, c(exports_common, "one_F"))
}

# =====================================================================
# Part S -- survival, vignette protocol
# =====================================================================
if ("S" %in% parts) {
  cat("\n=== Part S: survival (vignette protocol), T = 4 and T = 8 ===\n")
  dS <- file.path(out_dir, "S"); dir.create(dS, showWarnings = FALSE)
  THR <- list(T4 = c(4, 8, 14, 20), T8 = c(3, 5, 7, 9, 11, 14, 17, 20))
  # The vignette's grid c(1e-5, 1e-4, 1e-3) is too small for any penalty to
  # bite: in the smoke run every arm (lasso included) agreed to +-0.0007
  # C-index. Extended upward for every arm alike, so the penalties can differ.
  S_LAMBDAS <- 10^(-5:-1)
  one_S <- function(task) {
    suppressMessages(library(orthoMTL))
    setting <- task$setting; seed <- task$seed
    f <- file.path(out_dir, "S", sprintf("%s__seed%03d.rds", setting, seed))
    if (file.exists(f)) return(invisible(f))
    thr <- THR[[setting]]; Tt <- length(thr)
    sim <- simulate_mtl(n = 600, p = 15, n_signals = 5, thresholds = thr,
                        effect_strength = 1.2, seed = 1000 * Tt + seed)
    Y <- create_longitudinal_labels(sim$SurvTime, sim$Event, thr)
    W <- create_indicator_matrix(Y)
    tr <- 1:300; te <- 301:600
    X <- sim$X; folds <- rep(1:5, length.out = length(tr))
    Kd <- create_constraint_matrix(Tt)
    gt <- sim$ground_truth
    signal <- as.numeric(colnames(X) %in% gt$signal_features)
    auc <- function(score, lab) {                       # Mann-Whitney AUC
      r <- rank(score); n1 <- sum(lab == 1); n0 <- sum(lab == 0)
      (sum(r[lab == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    }
    # vignette grid; arms differ only in penalty / K / alpha set
    arms <- list(
      V_K   = list(penalty = "K",       K = Kd,           alphas = c(0, 0.5, 1), diag = c(0.5, Tt, 2 * Tt)),
      V_O   = list(penalty = "omega_O", K = NULL,         alphas = c(0, 0.5, 1), diag = 0.5),
      K0    = list(penalty = "K",       K = Kd,           alphas = 0,            diag = c(0.5, Tt, 2 * Tt)),
      O0    = list(penalty = "omega_O", K = NULL,         alphas = 0,            diag = 0.5),
      I0    = list(penalty = "K",       K = diag(1, Tt),  alphas = 0,            diag = 1),
      L     = list(penalty = "K",       K = diag(1, Tt),  alphas = 1,            diag = 1))
    rows <- lapply(names(arms), function(a) {
      arm <- arms[[a]]
      t0 <- Sys.time()
      cv <- cv_orthoMTL(X[tr, ], Y[tr, ], W[tr, ], K = arm$K,
                        lambdas = S_LAMBDAS, alphas = arm$alphas,
                        stepsizes = c(1, 2, 5), diag_vals = arm$diag,
                        survival = TRUE, penalty = arm$penalty, folds = folds,
                        n_cores = 1, seed = 42, verbose = FALSE)
      b <- cv$best
      Kf <- if (is.null(arm$K)) NULL else { k <- arm$K; diag(k) <- b$diag_val; k }
      fit <- orthoMTL(X[tr, ], Y[tr, ], lambda = b$lambda, alpha = b$alpha,
                      step_size = b$stepsize, K = Kf, penalty = arm$penalty,
                      survival = TRUE, censored.mat = W[tr, ], seed = 42)
      pred <- t(apply(X[te, ] %*% fit$B, 1, nnmaxheap_C))
      data.frame(setting = setting, seed = seed, strategy = a,
                 lambda = b$lambda, alpha = b$alpha, diag_val = b$diag_val,
                 cv = b$cv_score, cindex = cindex_mtl(Y[te, ], pred),
                 support_auc = auc(sqrt(rowSums(fit$B^2)), signal),
                 secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
    })
    save_atomic(do.call(rbind, rows), f)
    invisible(f)
  }
  grid <- expand.grid(setting = names(THR), seed = seq_len(n_seeds), stringsAsFactors = FALSE)
  tasks <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))
  pend <- Filter(function(t) !file.exists(file.path(dS, sprintf("%s__seed%03d.rds",
                   t$setting, t$seed))), tasks)
  cat(sprintf("  %d tasks, %d pending\n", length(tasks), length(pend)))
  run_tasks(pend, one_S, c(exports_common, "one_S", "THR", "S_LAMBDAS"))
}

# =====================================================================
# Aggregation and report
# =====================================================================
report_reg <- function(res, scen_list, title, strategies) {
  agg <- paired(res)
  cat(sprintf("\n--- %s: paired Delta test-MSE vs A (identity, lambda by CV); < 0 = better ---\n", title))
  cat("harm = share of seeds worse than A; worst = largest single-seed loss\n")
  for (sc in names(scen_list)) for (nz in scen_list[[sc]]$noise) {
    a <- agg[agg$scenario == sc & agg$noise == nz, ]
    if (nrow(a) == 0L) next
    cat(sprintf("%-12s noise %-4g  A MSE = %.3f  (n = %d)\n", sc, nz, a$base[1], a$n[1]))
    for (st in strategies) {
      z <- a[a$strategy == st, ]; if (nrow(z) == 0L) next
      cat(sprintf("   %-12s Delta=%+.4f +- %.4f  harm=%3.0f%%  worst=%+.4f  (%s)\n",
                  st, z$delta, z$ci, 100 * z$harm, z$worst, verdict(z)))
    }
  }
  # head-to-head O vs C (best K practice, CV over K) and O vs L
  pairs <- list(c("O", "C"), c("O", "L"))
  if ("OG" %in% res$strategy) pairs <- c(pairs, list(c("OG", "C"), c("OG", "O"), c("OG", "L")))
  for (pair in pairs) {
    r2 <- res[res$strategy %in% pair, ]
    r2$strategy <- ifelse(r2$strategy == pair[2], "A", "X")   # reuse paired(): base = pair[2]
    h <- paired(r2)
    cat(sprintf("\n  head-to-head %s - %s (< 0: %s better):\n", pair[1], pair[2], pair[1]))
    for (i in seq_len(nrow(h)))
      cat(sprintf("   %-12s noise %-4g Delta=%+.4f +- %.4f  (%s)\n", h$scenario[i], h$noise[i],
                  h$delta[i], h$ci[i], verdict(h[i, ])))
  }
  agg
}

# =====================================================================
# Part OG -- omega_O with the gamma ridge knob (penalty = "omega_O", gamma),
# (gamma, lambda) by the same 3-fold inner CV, on the R and F tasks. Data are
# regenerated with the exact streams of parts R and F; each task first refits
# the stored arm A and stops on any mismatch (pairing check). Needs the gamma
# build of orthoMTL. Results in ABO_OUT/OG_R and ABO_OUT/OG_F, merged into the
# R and F reports below as strategy "OG".
# =====================================================================
GAMMAS <- c(0, 0.25, 1, 4, 16)
if ("OG" %in% parts) {
  cat("\n=== Part OG: omega_O + gamma ridge on the R and F tasks ===\n")
  stopifnot("Part OG needs the gamma build (gamma= argument)" =
              "gamma" %in% names(formals(orthoMTL)))
  dir.create(file.path(out_dir, "OG_R"), showWarnings = FALSE)
  dir.create(file.path(out_dir, "OG_F"), showWarnings = FALSE)
  one_OG <- function(task) {
    suppressMessages(library(orthoMTL))
    part <- task$part; scen <- task$scenario; noise <- task$noise; seed <- task$seed
    f <- file.path(out_dir, paste0("OG_", part), sprintf("%s__noise%s__seed%03d.rds", scen, noise, seed))
    if (file.exists(f)) return(invisible(f))
    if (part == "R") {                                 # stream of part R (= K-01c)
      cfg <- SCEN_R[[scen]]
      set.seed(1e6 * match(scen, names(SCEN_R)) + 1e4 * match(noise, cfg$noise) + seed)
      Btrue <- make_B(scen, cfg)
      stored <- readRDS(file.path(k01c_dir, sprintf("%s__noise%s__seed%03d.rds", scen, noise, seed)))
    } else {                                           # stream of part F
      cfg <- SCEN_F[[scen]]
      set.seed(2e6 + 1e5 * match(scen, names(SCEN_F)) + 1e4 * match(noise, cfg$noise) + seed)
      Btrue <- make_B_F(scen, cfg)
      stored <- readRDS(file.path(out_dir, "F", sprintf("%s__noise%s__seed%03d.rds", scen, noise, seed)))
    }
    Xtr <- matrix(rnorm(cfg$n * cfg$p), cfg$n, cfg$p)
    Xte <- matrix(rnorm(NTEST * cfg$p), NTEST, cfg$p)
    Ytr <- Xtr %*% Btrue + matrix(rnorm(cfg$n * cfg$T), cfg$n) * noise
    Yte <- Xte %*% Btrue
    folds <- sample(rep(seq_len(CV_K), length.out = cfg$n))
    a <- stored[stored$strategy == "A", ]
    mseA <- fit_mse(Xtr, Ytr, Xte, Yte, a$lambda, list(penalty = "K", K = diag(1, cfg$T), alpha = 0))
    if (!isTRUE(all.equal(mseA, a$mse, tolerance = 1e-10)))
      stop(sprintf("OG pairing check failed for %s: refit A %.12g vs stored %.12g", basename(f), mseA, a$mse))
    cv <- outer(GAMMAS, LAMBDAS, Vectorize(function(g, l) mean(vapply(seq_len(CV_K), function(k) {
      tr <- folds != k; va <- folds == k
      fit_mse(Xtr[tr, , drop = FALSE], Ytr[tr, , drop = FALSE], Xtr[va, , drop = FALSE],
              Ytr[va, , drop = FALSE], l, list(penalty = "omega_O", K = NULL, alpha = 0, gamma = g))
    }, numeric(1)))))
    stopifnot("an OG fit failed (NA cv): check the library" = all(is.finite(cv)))
    ij <- arrayInd(which.min(cv), dim(cv))
    g <- GAMMAS[ij[1]]; l <- LAMBDAS[ij[2]]
    mse <- fit_mse(Xtr, Ytr, Xte, Yte, l, list(penalty = "omega_O", K = NULL, alpha = 0, gamma = g))
    save_atomic(data.frame(scenario = scen, noise = noise, seed = seed, strategy = "OG",
                           gamma = g, lambda = l, cv = cv[ij], mse = mse,
                           cv_gamma0 = min(cv[1, ])), f)
    invisible(f)
  }
  grid <- rbind(
    do.call(rbind, lapply(names(SCEN_R), function(s) expand.grid(part = "R", scenario = s,
      noise = SCEN_R[[s]]$noise, seed = seq_len(n_seeds), stringsAsFactors = FALSE))),
    do.call(rbind, lapply(names(SCEN_F), function(s) expand.grid(part = "F", scenario = s,
      noise = SCEN_F[[s]]$noise, seed = seq_len(n_seeds), stringsAsFactors = FALSE))))
  grid <- grid[order(grid$seed), ]
  tasks <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))
  pend <- Filter(function(t) !file.exists(file.path(out_dir, paste0("OG_", t$part),
                   sprintf("%s__noise%s__seed%03d.rds", t$scenario, t$noise, t$seed))), tasks)
  cat(sprintf("  %d tasks, %d pending\n", length(tasks), length(pend)))
  run_tasks(pend, one_OG, c(exports_common, "one_OG", "GAMMAS"))
}
og_rows <- function(part) {
  fs <- list.files(file.path(out_dir, paste0("OG_", part)), pattern = "\\.rds$", full.names = TRUE)
  if (length(fs) == 0L) return(NULL)
  do.call(rbind, lapply(fs, readRDS))
}

out <- list()
if ("R" %in% parts) {
  fs <- list.files(file.path(out_dir, "R"), pattern = "\\.rds$", full.names = TRUE)
  if (length(fs)) {
    res <- do.call(rbind, lapply(fs, readRDS))
    og <- og_rows("R")
    if (!is.null(og)) res <- rbind(res, og[, names(res)])
    cat(sprintf("\nPart R: %d of %d tasks finished (pairing check passed in each)\n",
                length(fs), 4 * 3 * n_seeds))
    out$R <- list(res = res, agg = report_reg(res, SCEN_R, "Part R (K-01c scenarios)",
                  c("B_uniform", "B_diffusion", "B_kernel", "C", "D", "E", "O", "OG", "L")))
  }
}
if ("F" %in% parts) {
  fs <- list.files(file.path(out_dir, "F"), pattern = "\\.rds$", full.names = TRUE)
  if (length(fs)) {
    res <- do.call(rbind, lapply(fs, readRDS))
    og <- og_rows("F")
    if (!is.null(og)) { og$shape <- NA; og$rho <- NA; res <- rbind(res, og[, names(res)]) }
    cat(sprintf("\nPart F: %d of %d tasks finished\n", length(fs), 3 * 3 * n_seeds))
    out$F <- list(res = res, agg = report_reg(res, SCEN_F, "Part F (new scenarios)",
                  c("B_uniform", "B_diffusion", "C", "O", "OG", "L")))
    cat("\n  CV-selected (shape, rho) under C:\n")
    sel <- res[res$strategy == "C", ]
    print(round(prop.table(table(sel$scenario, paste(sel$shape, ifelse(is.na(sel$rho), "", sel$rho))), 1), 2))
  }
}
if ("S" %in% parts) {
  fs <- list.files(file.path(out_dir, "S"), pattern = "\\.rds$", full.names = TRUE)
  if (length(fs)) {
    res <- do.call(rbind, lapply(fs, readRDS))
    res$scenario <- res$setting; res$noise <- 0
    cat(sprintf("\nPart S: %d of %d tasks finished\n", length(fs), 2 * n_seeds))
    for (m in c("cindex", "support_auc")) {
      agg <- paired(res, base = "V_K", keys = c("scenario", "noise", "seed"), value = m)
      cat(sprintf("\n--- Part S: paired Delta %s vs V_K (vignette protocol with diffusion K); > 0 = better ---\n", m))
      for (s in unique(agg$scenario)) {
        a <- agg[agg$scenario == s, ]
        cat(sprintf("%s  V_K mean %s = %.4f (n = %d)\n", s, m, a$base[1], a$n[1]))
        for (i in seq_len(nrow(a)))
          cat(sprintf("   %-6s Delta=%+.4f +- %.4f  (%s)\n", a$strategy[i], a$delta[i], a$ci[i],
                      verdict(a[i, ], lower_is_better = FALSE)))
      }
      out$S[[m]] <- agg
    }
    cat("\n  CV-selected alpha under V_O / V_K (share):\n")
    print(round(prop.table(table(res$strategy[res$strategy %in% c("V_K", "V_O")],
                                 res$alpha[res$strategy %in% c("V_K", "V_O")]), 1), 2))
    cat(sprintf("\n  mean seconds per arm (CV + refit): %s\n",
                paste(sprintf("%s %.0f", names(tapply(res$secs, res$strategy, mean)),
                              tapply(res$secs, res$strategy, mean)), collapse = ", ")))
    out$S$res <- res
  }
}
for (part in c("R", "F")) {
  og <- og_rows(part)
  if (is.null(og)) next
  cat(sprintf("\n  Part OG on %s: CV-selected gamma (share of tasks)\n", part))
  print(round(prop.table(table(og$scenario, og$gamma), 1), 2))
}
if (length(out)) saveRDS(out, file.path(out_dir, "omega_o_aggregate.rds"))
cat(sprintf("\nResults: %s\n", out_dir))
