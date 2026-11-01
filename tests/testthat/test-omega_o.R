# Tests for the omega_O penalty: the max-cut SDP helper and its use in orthoMTL()

test_that(".maxcut_sdp matches the closed form for two tasks", {
  # T = 2: max <G, Z> over correlation matrices = G11 + G22 + 2|G12|
  W <- matrix(c(1, 2, 0, -1, 1, 3), ncol = 2)
  G <- crossprod(W)
  res <- orthoMTL:::.maxcut_sdp(G)
  expect_equal(res$value, G[1, 1] + G[2, 2] + 2 * abs(G[1, 2]), tolerance = 1e-9)
  expect_lt(res$gap, 1e-8)
})

test_that(".maxcut_sdp equals sum(G) on balanced Gram matrices (Proposition C)", {
  set.seed(1)
  W <- abs(matrix(rnorm(8 * 5), 8, 5))       # all inner products >= 0
  flip <- c(1, -1, 1, -1, -1)
  W <- W %*% diag(flip)                     # balanced, not entrywise >= 0
  G <- crossprod(W)
  res <- orthoMTL:::.maxcut_sdp(G)
  expect_equal(res$value, sum(abs(G)), tolerance = 1e-8)
  expect_lt(res$gap, 1e-6)
})

test_that(".maxcut_sdp is bracketed by max-cut and the pairwise box bound", {
  set.seed(2)
  for (rep in 1:5) {
    W <- matrix(rnorm(6 * 5), 6, 5)
    G <- crossprod(W)
    res <- orthoMTL:::.maxcut_sdp(G)
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), 5)))
    maxcut <- max(apply(signs, 1, function(e) sum(G * tcrossprod(e))))
    expect_gte(res$value, maxcut - 1e-8)            # elliptope contains every eps eps^T
    expect_lte(res$value, sum(abs(G)) + 1e-8)       # elliptope inside the box |Z_ij| <= 1
    expect_lt(res$gap, 1e-5 * max(1, res$value))    # certified near-optimal
    expect_equal(diag(res$Z), rep(1, 5), tolerance = 1e-12)
    expect_gte(min(eigen(res$Z, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
  }
})

test_that(".maxcut_sdp gives the frustrated 120-degree value", {
  # three unit vectors at 120 degrees: pairwise box bound 6, SDP 4.5, max-cut 4
  ang <- c(0, 2 * pi / 3, 4 * pi / 3)
  W <- rbind(cos(ang), sin(ang))
  res <- orthoMTL:::.maxcut_sdp(crossprod(W))
  expect_equal(res$value, 4.5, tolerance = 1e-8)
  # a collapsed rank-1 warm start (optimal for a balanced G) must not trap it at max-cut (4)
  V1 <- matrix(c(1, 0, 0), 3, 3)
  res1 <- orthoMTL:::.maxcut_sdp(crossprod(W), V0 = V1)
  expect_equal(res1$value, 4.5, tolerance = 1e-6)
})

test_that("penalty = 'K' (default) is unchanged by the new argument", {
  set.seed(42)
  X <- matrix(rnorm(60 * 6), 60, 6)
  Y <- X %*% matrix(rnorm(6 * 3), 6, 3) + matrix(rnorm(60 * 3), 60)
  K <- create_constraint_matrix(3, diag_val = 2)
  f0 <- orthoMTL(X, Y, lambda = 1e-2, K = K, max_iter = 2000)
  f1 <- orthoMTL(X, Y, lambda = 1e-2, K = K, max_iter = 2000, penalty = "K")
  expect_identical(f0$B, f1$B)
  expect_identical(f0$obj, f1$obj)
  expect_identical(f0$hyperparameters$penalty, "K")
  expect_true(is.na(f0$sdp_gap))
})

test_that("penalty = 'omega_O' fits, reports a certified objective, ignores K", {
  set.seed(42)
  X <- matrix(rnorm(80 * 6), 80, 6)
  Y <- X %*% qr.Q(qr(matrix(rnorm(6 * 3), 6, 3))) + 0.1 * matrix(rnorm(80 * 3), 80)
  f <- orthoMTL(X, Y, lambda = 1e-2, penalty = "omega_O", schedule = "log", max_iter = 3000)
  expect_s3_class(f, "orthoMTL")
  expect_identical(f$hyperparameters$penalty, "omega_O")
  expect_lt(f$sdp_gap, 1e-5)
  # the reported objective uses the SDP value at B
  pen <- orthoMTL:::.maxcut_sdp(crossprod(f$B))$value
  obj <- 0.5 * sum((X %*% f$B - Y)^2) / nrow(X) + 0.5 * 1e-2 * pen
  expect_equal(f$obj, obj, tolerance = 1e-6)
  fK <- orthoMTL(X, Y, lambda = 1e-2, penalty = "omega_O", schedule = "log", max_iter = 3000,
                 K = create_constraint_matrix(3))
  expect_identical(f$B, fK$B)
  expect_output(print(f), "omega_O \\(SDP gap")
})

test_that("penalty = 'omega_O' guards: disjoint refused, p < T warns", {
  X <- matrix(rnorm(30 * 2), 30, 2)
  Y <- matrix(rnorm(30 * 3), 30, 3)
  expect_error(orthoMTL(X, Y, penalty = "omega_O", disjoint = TRUE), "disjoint")
  expect_warning(orthoMTL(X, Y, lambda = 1e-2, penalty = "omega_O", max_iter = 50),
                 "fewer features than tasks")
})

test_that("gamma: validated, a no-op at 0, enters the objective, tends to ridge", {
  set.seed(5)
  X <- matrix(rnorm(80 * 6), 80, 6)
  Y <- X %*% matrix(rnorm(6 * 3), 6, 3) + matrix(rnorm(80 * 3), 80)
  expect_error(orthoMTL(X, Y, gamma = 1), "'omega_O' only")
  expect_error(orthoMTL(X, Y, penalty = "omega_O", gamma = -1), "non-negative")
  f0 <- orthoMTL(X, Y, lambda = 1e-2, penalty = "omega_O", schedule = "log", max_iter = 2000)
  fz <- orthoMTL(X, Y, lambda = 1e-2, penalty = "omega_O", gamma = 0, schedule = "log", max_iter = 2000)
  expect_identical(f0$B, fz$B)
  fg <- orthoMTL(X, Y, lambda = 1e-2, penalty = "omega_O", gamma = 2, schedule = "log", max_iter = 3000)
  expect_identical(fg$hyperparameters$gamma, 2)
  pen <- orthoMTL:::.maxcut_sdp(crossprod(fg$B))$value + 2 * sum(fg$B^2)
  expect_equal(fg$obj, 0.5 * sum((X %*% fg$B - Y)^2) / nrow(X) + 0.5 * 1e-2 * pen,
               tolerance = 1e-6)
  expect_output(print(fg), "gamma ridge 2")
  # gamma -> infinity with lambda * gamma fixed: Omega_O^2 <= T ||W||_F^2 becomes
  # negligible next to gamma ||W||_F^2, so the fit tends to ridge (K = I)
  fr <- orthoMTL(X, Y, lambda = 0.5, K = diag(3), schedule = "log", max_iter = 5000)
  fbig <- orthoMTL(X, Y, lambda = 0.5 / 1e4, penalty = "omega_O", gamma = 1e4,
                   schedule = "log", max_iter = 5000)
  expect_lt(max(abs(fbig$B - fr$B)) / max(abs(fr$B)), 2e-2)
})

test_that("cv_orthoMTL searches gammas for omega_O only", {
  set.seed(4)
  X <- matrix(rnorm(40 * 5), 40, 5)
  Y <- X %*% matrix(rnorm(5 * 3), 5, 3) + matrix(rnorm(40 * 3), 40)
  cvO <- cv_orthoMTL(X, Y, lambdas = c(1e-3, 1e-2), stepsizes = 0.1, survival = FALSE,
                     penalty = "omega_O", gammas = c(0, 1), folds = rep(1:2, 20),
                     n_cores = 1, verbose = FALSE)
  expect_equal(cvO$n_configs, 4L)
  expect_true(cvO$best$gamma %in% c(0, 1))
  expect_true("gamma" %in% names(cvO$results))
  cvK <- cv_orthoMTL(X, Y, lambdas = 1e-2, stepsizes = 0.1, diag_vals = 1, survival = FALSE,
                     gammas = c(0, 1), folds = rep(1:2, 20), n_cores = 1, verbose = FALSE)
  expect_null(cvK$best$gamma)
  expect_false("gamma" %in% names(cvK$results))
  expect_error(cv_orthoMTL(X, Y, survival = FALSE, penalty = "omega_O", gammas = -1,
                           folds = rep(1:2, 20), n_cores = 1, verbose = FALSE),
               "non-negative")
})

test_that("cv_orthoMTL threads penalty and collapses diag_vals for omega_O", {
  set.seed(3)
  X <- matrix(rnorm(40 * 5), 40, 5)
  Y <- X %*% matrix(rnorm(5 * 3), 5, 3) + matrix(rnorm(40 * 3), 40)
  cv <- cv_orthoMTL(X, Y, lambdas = c(1e-3, 1e-2), stepsizes = 0.1,
                    diag_vals = c(0.5, 1, 2), survival = FALSE,
                    penalty = "omega_O", folds = rep(1:2, 20),
                    n_cores = 1, verbose = FALSE)
  expect_equal(cv$n_configs, 2L)
})
