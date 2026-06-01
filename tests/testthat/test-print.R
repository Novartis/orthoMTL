# Tests for print, summary, and coef S3 methods

test_that("print.orthoMTL produces output and returns invisibly", {
  set.seed(42)
  X <- matrix(rnorm(60), 12, 5)
  Y <- matrix(rnorm(36), 12, 3)
  fit <- orthoMTL(X, Y, lambda = 0.01)

  expect_output(print(fit), "orthoMTL model")
  expect_output(print(fit), "Converged")
  expect_output(print(fit), "lambda")

  out <- print(fit)
  expect_s3_class(out, "orthoMTL")
})

test_that("print.orthoMTL shows correct mode", {
  X <- matrix(rnorm(60), 12, 5)

  # Regression
  Y_reg <- matrix(rnorm(36), 12, 3)
  fit_reg <- orthoMTL(X, Y_reg, lambda = 0.01)
  expect_output(print(fit_reg), "regression")

  # Classification
  Y_cls <- matrix(sample(c(-1, 1), 36, replace = TRUE), 12, 3)
  fit_cls <- orthoMTL(X, Y_cls, lambda = 0.01, logistic = TRUE)
  expect_output(print(fit_cls), "classification")
})

test_that("summary.orthoMTL produces detailed output", {
  set.seed(42)
  X <- matrix(rnorm(60), 12, 5)
  colnames(X) <- paste0("V", 1:5)
  Y <- matrix(rnorm(36), 12, 3)
  fit <- orthoMTL(X, Y, lambda = 0.01)

  expect_output(summary(fit), "orthoMTL model")
  expect_output(summary(fit), "Hyperparameters")
  expect_output(summary(fit), "Coefficients")
  expect_output(summary(fit), "Sparsity")
  expect_output(summary(fit), "Top")
})

test_that("summary.orthoMTL shows task thresholds for survival", {
  set.seed(42)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  SurvTime <- rexp(n, 0.1)
  Event <- rbinom(n, 1, 0.7)
  Y <- create_longitudinal_labels(SurvTime, Event, c(4, 6, 10))
  W <- create_indicator_matrix(Y)

  fit <- orthoMTL(X, Y, lambda = 0.01, survival = TRUE, censored.mat = W)

  expect_output(summary(fit), "Task thresholds")
})

test_that("coef.orthoMTL returns the coefficient matrix", {
  set.seed(42)
  X <- matrix(rnorm(60), 12, 5)
  Y <- matrix(rnorm(36), 12, 3)
  fit <- orthoMTL(X, Y, lambda = 0.01)

  B <- coef(fit)
  expect_true(is.matrix(B))
  expect_identical(B, fit$B)
})
