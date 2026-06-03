# Smoke tests for cv_orthoMTL

test_that("cv_orthoMTL returns S3 class with expected structure", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  SurvTime <- rexp(n, 0.1)
  Event <- rbinom(n, 1, 0.7)

  Y <- create_longitudinal_labels(SurvTime, Event, c(4, 6, 10))
  W <- create_indicator_matrix(Y)
  K <- create_constraint_matrix(n_tasks)
  folds <- rep(1:2, length.out = n)

  cv_res <- cv_orthoMTL(
    X.train = X, Y.train = Y, W.train = W, K = K,
    lambdas = c(1e-3), alphas = 0,
    stepsizes = c(0.1), diag_vals = c(0.5),
    survival = TRUE, disjoint = FALSE,
    folds = folds, n_cores = 1, seed = 42, verbose = FALSE
  )

  expect_s3_class(cv_res, "cv_orthoMTL")
  expect_true(is.list(cv_res$best))
  expect_true(is.data.frame(cv_res$results))
  expect_equal(length(cv_res$folds), n)
  expect_true("cv_score" %in% colnames(cv_res$results))
  expect_equal(cv_res$n_folds, 2)
})

test_that("cv_orthoMTL best has all expected fields", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  folds <- rep(1:2, length.out = n)

  cv_res <- cv_orthoMTL(
    X.train = X, Y.train = Y,
    lambdas = c(0.01), stepsizes = c(0.1), diag_vals = c(1),
    survival = FALSE, folds = folds, n_cores = 1, verbose = FALSE
  )

  expect_true(all(c("lambda", "alpha", "stepsize", "diag_val", "cv_score")
                  %in% names(cv_res$best)))
})

test_that("cv_orthoMTL searches multiple configurations", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  folds <- rep(1:2, length.out = n)

  cv_res <- cv_orthoMTL(
    X.train = X, Y.train = Y,
    lambdas = c(0.001, 0.01), stepsizes = c(0.1, 0.5),
    diag_vals = c(0.5, 1),
    survival = FALSE, folds = folds, n_cores = 1, verbose = FALSE
  )

  # 2 lambdas x 2 stepsizes x 2 diag_vals x 1 alpha = 8
  expect_equal(cv_res$n_configs, 8)
  expect_equal(nrow(cv_res$results), 8)

  # Results should be sorted by cv_score descending
  expect_true(all(diff(cv_res$results$cv_score) <= 0))
})

test_that("cv_orthoMTL generates folds with warning when not provided", {
  set.seed(42)
  n <- 20; p <- 3; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  expect_warning(
    cv_orthoMTL(X.train = X, Y.train = Y,
                lambdas = c(0.01), stepsizes = c(0.1), diag_vals = c(1),
                survival = FALSE, folds = NULL, n_cores = 1, verbose = FALSE),
    "not provided"
  )
})

test_that("cv_orthoMTL validates inputs", {
  expect_error(
    cv_orthoMTL(X.train = NULL, Y.train = NULL,
                lambdas = 1, stepsizes = 1, diag_vals = 1),
    "must be provided"
  )
})

test_that("print.cv_orthoMTL produces output", {
  set.seed(42)
  n <- 20; p <- 3; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  folds <- rep(1:2, length.out = n)

  cv_res <- cv_orthoMTL(
    X.train = X, Y.train = Y,
    lambdas = c(0.01), stepsizes = c(0.1), diag_vals = c(1),
    survival = FALSE, folds = folds, n_cores = 1, verbose = FALSE
  )

  expect_output(print(cv_res), "Cross-validation")
  expect_output(print(cv_res), "Best configuration")
})
