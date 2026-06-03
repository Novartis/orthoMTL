# Tests for the core orthoMTL solver

test_that("orthoMTL returns S3 class with expected structure", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  colnames(Y) <- paste0("T", seq_len(n_tasks))

  fit <- orthoMTL(X, Y, lambda = 0.01)

  expect_s3_class(fit, "orthoMTL")
  expect_true(is.matrix(fit$B))
  expect_equal(dim(fit$B), c(p, n_tasks))
  expect_identical(fit$B, fit$W)  # alias
  expect_true(is.numeric(fit$obj))
  expect_true(is.numeric(fit$imax))
  expect_true(is.logical(fit$converged))
  expect_true(is.list(fit$hyperparameters))
  expect_equal(fit$n_features, p)
  expect_equal(fit$n_tasks, n_tasks)
  expect_equal(fit$n_obs, n)
  expect_equal(fit$feature_names, colnames(X))
  expect_equal(fit$task_names, colnames(Y))
})

test_that("orthoMTL preserves row/column names on B", {
  X <- matrix(rnorm(40), 8, 5)
  colnames(X) <- paste0("f", 1:5)
  Y <- matrix(rnorm(24), 8, 3)
  colnames(Y) <- c("t1", "t2", "t3")

  fit <- orthoMTL(X, Y, lambda = 0.01)

  expect_equal(rownames(fit$B), colnames(X))
  expect_equal(colnames(fit$B), colnames(Y))
})

test_that("orthoMTL converges on simple regression problem", {
  set.seed(42)
  n <- 100; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  W_true <- matrix(rnorm(p * n_tasks), p, n_tasks)
  Y <- X %*% W_true + matrix(rnorm(n * n_tasks), n) * 0.01

  fit <- orthoMTL(X, Y, lambda = 1e-4, step_size = 0.1)

  expect_true(fit$converged)
  expect_true(fit$imax < 1e6)
})

test_that("orthoMTL respects logistic mode", {
  set.seed(42)
  n <- 50; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(sample(c(-1, 1), n * n_tasks, replace = TRUE), n, n_tasks)

  fit <- orthoMTL(X, Y, lambda = 0.01, logistic = TRUE)

  expect_s3_class(fit, "orthoMTL")
  expect_true(fit$hyperparameters$logistic)
})

test_that("logistic gradient is correct under label flip (REG-01)", {
  # The logistic loss log(1+exp(-Y*Xw)) is invariant under (Y, W) -> (-Y, -W),
  # so flipping all labels must negate the fitted coefficients. The old
  # gradient (-(Y - sigmoid(Xw)), a y in {0,1} form) breaks this symmetry.
  set.seed(7)
  n <- 60; p <- 4; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(sample(c(-1, 1), n * n_tasks, replace = TRUE), n, n_tasks)

  fit_pos <- orthoMTL(X,  Y, lambda = 0.01, logistic = TRUE)
  fit_neg <- orthoMTL(X, -Y, lambda = 0.01, logistic = TRUE)

  expect_equal(fit_neg$B, -fit_pos$B, tolerance = 1e-6)
})

test_that("logistic loss does not overflow on large scores (REG-02)", {
  # Large |Xw| made the naive log(1+exp(-Y*Xw)) overflow to Inf/NaN.
  set.seed(11)
  n <- 40; p <- 3; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p) * 100
  Y <- matrix(sample(c(-1, 1), n * n_tasks, replace = TRUE), n, n_tasks)

  fit <- orthoMTL(X, Y, lambda = 0.01, logistic = TRUE, step_size = 1)

  expect_true(is.finite(fit$obj))
  expect_false(anyNA(fit$B))
})

test_that("orthoMTL survival mode requires censored.mat", {
  X <- matrix(rnorm(20), 4, 5)
  Y <- matrix(c(1, 0, NA, 1, 1, 1, NA, 0), 4, 2)

  expect_error(
    orthoMTL(X, Y, survival = TRUE, censored.mat = NULL),
    "censoring information"
  )
})

test_that("orthoMTL survival mode runs with valid censored.mat", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  SurvTime <- rexp(n, 0.1)
  Event <- rbinom(n, 1, 0.7)

  Y <- create_longitudinal_labels(SurvTime, Event, c(4, 6, 10))
  W <- create_indicator_matrix(Y)

  fit <- orthoMTL(X, Y, lambda = 0.01, survival = TRUE, censored.mat = W)

  expect_s3_class(fit, "orthoMTL")
  expect_true(fit$hyperparameters$survival)
})

test_that("orthoMTL disjoint mode runs", {
  set.seed(42)
  n <- 50; p <- 5; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  fit <- orthoMTL(X, Y, lambda = 0.01, disjoint = TRUE, step_size = 5)

  expect_s3_class(fit, "orthoMTL")
  expect_true(fit$hyperparameters$disjoint)
})

test_that("orthoMTL elastic-net (alpha > 0) runs and stores alpha", {
  set.seed(42)
  n <- 50; p <- 5; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  fit <- orthoMTL(X, Y, lambda = 0.01, alpha = 0.5)

  expect_s3_class(fit, "orthoMTL")
  expect_equal(fit$hyperparameters$alpha, 0.5)
})

test_that("orthoMTL rejects alpha outside [0, 1]", {
  X <- matrix(rnorm(20), 4, 5)
  Y <- matrix(rnorm(8), 4, 2)
  expect_error(orthoMTL(X, Y, alpha = 1.5), "\\[0, 1\\]")
  expect_error(orthoMTL(X, Y, alpha = -0.1), "\\[0, 1\\]")
})

test_that("elastic-net mixing endpoints both produce finite fits", {
  set.seed(1)
  n <- 40; p <- 6; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  fit0 <- orthoMTL(X, Y, lambda = 0.01, alpha = 0)   # pure orthogonality
  fit1 <- orthoMTL(X, Y, lambda = 0.01, alpha = 1)   # pure Lasso

  expect_false(anyNA(fit0$B))
  expect_false(anyNA(fit1$B))
  expect_true(is.finite(fit0$obj))
  expect_true(is.finite(fit1$obj))
})

test_that("orthoMTL stops on invalid inputs", {
  expect_error(orthoMTL(NULL, NULL), "need to be provided")
  expect_error(orthoMTL(matrix(1, 2, 2), matrix(1, 2, 2), lambda = -1),
               "positive")
})

test_that("orthoMTL tol parameter is used", {
  set.seed(42)
  n <- 50; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  fit_strict <- orthoMTL(X, Y, lambda = 0.01, tol = 1e-8)
  fit_loose <- orthoMTL(X, Y, lambda = 0.01, tol = 1e-2)

  # Stricter tolerance should take more (or equal) iterations
  expect_true(fit_strict$imax >= fit_loose$imax)
})

test_that("orthoMTL W_0 parameter works", {
  set.seed(42)
  n <- 50; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  W_0 <- matrix(0.1, p, n_tasks)
  fit <- orthoMTL(X, Y, lambda = 0.01, W_0 = W_0)

  expect_s3_class(fit, "orthoMTL")
})
