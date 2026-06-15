# Tests for mode-aware CV scoring (GEN-02)

test_that("cv_orthoMTL defaults to rmse for plain regression", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  folds <- rep(1:2, length.out = n)

  cv_res <- cv_orthoMTL(
    X.train = X, Y.train = Y,
    lambdas = c(0.001, 0.01), stepsizes = c(0.1), diag_vals = c(1),
    survival = FALSE, folds = folds, n_cores = 1, verbose = FALSE
  )

  expect_equal(cv_res$metric, "rmse")
  # rmse: lower is better, so results are sorted ascending
  expect_true(all(diff(cv_res$results$cv_score) >= 0))
  expect_equal(cv_res$best$cv_score, min(cv_res$results$cv_score))
})

test_that("cv_orthoMTL defaults to auc for logistic mode", {
  set.seed(42)
  n <- 60; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  lin <- X[, 1] - X[, 2]
  Y <- matrix(0, n, n_tasks)
  for (k in seq_len(n_tasks)) {
    Y[, k] <- ifelse(stats::runif(n) < 1 / (1 + exp(-lin)), 1, -1)
  }
  folds <- rep(1:2, length.out = n)

  cv_res <- cv_orthoMTL(
    X.train = X, Y.train = Y,
    lambdas = c(0.001, 0.01), stepsizes = c(0.1), diag_vals = c(1),
    survival = FALSE, logistic = TRUE,
    folds = folds, n_cores = 1, verbose = FALSE
  )

  expect_equal(cv_res$metric, "auc")
  # auc: higher is better, so results are sorted descending
  expect_true(all(diff(cv_res$results$cv_score) <= 0))
  expect_equal(cv_res$best$cv_score, max(cv_res$results$cv_score))
})

test_that("cv_orthoMTL survival default stays cindex", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  SurvTime <- rexp(n, 0.1); Event <- rbinom(n, 1, 0.7)
  Y <- create_longitudinal_labels(SurvTime, Event, c(4, 6, 10))
  W <- create_indicator_matrix(Y)
  folds <- rep(1:2, length.out = n)

  cv_res <- cv_orthoMTL(
    X.train = X, Y.train = Y, W.train = W,
    lambdas = c(0.01), stepsizes = c(0.1), diag_vals = c(1),
    survival = TRUE, folds = folds, n_cores = 1, verbose = FALSE
  )
  expect_equal(cv_res$metric, "cindex")
})

test_that("cv_orthoMTL honours an explicit metric argument", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  folds <- rep(1:2, length.out = n)

  cv_res <- cv_orthoMTL(
    X.train = X, Y.train = Y, metric = "r2",
    lambdas = c(0.001, 0.01), stepsizes = c(0.1), diag_vals = c(1),
    survival = FALSE, folds = folds, n_cores = 1, verbose = FALSE
  )
  expect_equal(cv_res$metric, "r2")
  expect_true(all(diff(cv_res$results$cv_score) <= 0))  # higher better
})

test_that("cv_orthoMTL rejects unknown metric", {
  set.seed(42)
  X <- matrix(rnorm(60), 30, 2)
  Y <- matrix(rnorm(60), 30, 2)
  folds <- rep(1:2, length.out = 30)
  expect_error(
    cv_orthoMTL(X.train = X, Y.train = Y, metric = "bogus",
                lambdas = 0.01, stepsizes = 0.1, diag_vals = 1,
                survival = FALSE, folds = folds, n_cores = 1, verbose = FALSE),
    "'arg'"
  )
})

test_that("print.cv_orthoMTL shows the metric name", {
  set.seed(42)
  X <- matrix(rnorm(60), 30, 2)
  Y <- matrix(rnorm(60), 30, 2)
  folds <- rep(1:2, length.out = 30)
  cv_res <- cv_orthoMTL(
    X.train = X, Y.train = Y,
    lambdas = 0.01, stepsizes = 0.1, diag_vals = 1,
    survival = FALSE, folds = folds, n_cores = 1, verbose = FALSE
  )
  expect_output(print(cv_res), "CV rmse")
})
