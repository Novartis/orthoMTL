# Smoke tests for bootstrap_orthoMTL

test_that("bootstrap_orthoMTL returns S3 class with expected structure", {
  set.seed(42)
  n <- 20; p <- 4; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  colnames(Y) <- paste0("T", 1:n_tasks)

  boot_res <- bootstrap_orthoMTL(
    X = X, Y = Y, lambda = 0.01, step_size = 0.1,
    n_repeats = 3, n_cores = 1, verbose = FALSE
  )

  expect_s3_class(boot_res, "bootstrap_orthoMTL")
  expect_true(is.data.frame(boot_res$results))
  expect_true(is.list(boot_res$coefficients_real))
  expect_true(is.list(boot_res$coefficients_null))
  expect_equal(length(boot_res$coefficients_real), 3)
  expect_equal(length(boot_res$coefficients_null), 3)
  expect_equal(length(boot_res$obj_real), 3)
  expect_equal(length(boot_res$obj_null), 3)
  expect_equal(boot_res$n_repeats, 3)
  expect_equal(boot_res$n_features, p)
  expect_equal(boot_res$n_tasks, n_tasks)
})

test_that("bootstrap_orthoMTL results data.frame has expected columns", {
  set.seed(42)
  n <- 20; p <- 3; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  colnames(Y) <- paste0("T", 1:n_tasks)

  boot_res <- bootstrap_orthoMTL(
    X = X, Y = Y, lambda = 0.01,
    n_repeats = 2, n_cores = 1, verbose = FALSE
  )

  expected_cols <- c("id", "time", "coeff", "group", "repeat_id")
  expect_true(all(expected_cols %in% colnames(boot_res$results)))

  # Should have both groups
  expect_true(all(c("real", "null") %in% unique(boot_res$results$group)))

  # All features should appear
  expect_equal(sort(unique(boot_res$results$id)),
               sort(colnames(X)))
})

test_that("bootstrap_orthoMTL results have correct row count", {
  set.seed(42)
  n <- 20; p <- 3; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  colnames(Y) <- paste0("T", 1:n_tasks)
  n_repeats <- 4

  boot_res <- bootstrap_orthoMTL(
    X = X, Y = Y, lambda = 0.01,
    n_repeats = n_repeats, n_cores = 1, verbose = FALSE
  )

  # Expected: n_features * n_tasks * n_repeats * 2 (real + null)
  expected_rows <- p * n_tasks * n_repeats * 2
  expect_equal(nrow(boot_res$results), expected_rows)
})

test_that("bootstrap_orthoMTL works in survival mode", {
  set.seed(42)
  n <- 25; p <- 4; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  SurvTime <- rexp(n, 0.1)
  Event <- rbinom(n, 1, 0.7)

  Y <- create_longitudinal_labels(SurvTime, Event, c(4, 6, 10))
  W <- create_indicator_matrix(Y)
  K <- create_constraint_matrix(n_tasks)

  boot_res <- bootstrap_orthoMTL(
    X = X, Y = Y, lambda = 0.01, step_size = 0.1,
    K = K, survival = TRUE, censored.mat = W,
    n_repeats = 2, n_cores = 1, verbose = FALSE
  )

  expect_s3_class(boot_res, "bootstrap_orthoMTL")
  expect_equal(boot_res$n_tasks, n_tasks)
})

test_that("bootstrap_orthoMTL obj_real and obj_null are separate", {
  set.seed(42)
  n <- 20; p <- 3; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  boot_res <- bootstrap_orthoMTL(
    X = X, Y = Y, lambda = 0.01,
    n_repeats = 3, n_cores = 1, verbose = FALSE
  )

  # Real and null objectives should generally differ
  # (extremely unlikely to be identical with different data)
  expect_false(identical(boot_res$obj_real, boot_res$obj_null))
})

test_that("bootstrap_orthoMTL validates inputs", {
  expect_error(
    bootstrap_orthoMTL(X = NULL, Y = NULL),
    "must be provided"
  )

  X <- matrix(rnorm(20), 4, 5)
  Y <- matrix(rnorm(12), 4, 3)
  expect_error(
    bootstrap_orthoMTL(X = X, Y = Y, survival = TRUE, censored.mat = NULL),
    "censored.mat"
  )
})

test_that("print.bootstrap_orthoMTL produces output", {
  set.seed(42)
  n <- 20; p <- 3; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  boot_res <- bootstrap_orthoMTL(
    X = X, Y = Y, lambda = 0.01,
    n_repeats = 2, n_cores = 1, verbose = FALSE
  )

  expect_output(print(boot_res), "Bootstrap inference")
  expect_output(print(boot_res), "Real models")
  expect_output(print(boot_res), "Null models")
})
