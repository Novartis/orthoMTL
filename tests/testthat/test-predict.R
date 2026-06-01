# Tests for predict.orthoMTL

test_that("predict.orthoMTL returns correct dimensions", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  colnames(Y) <- paste0("T", 1:n_tasks)

  fit <- orthoMTL(X, Y, lambda = 0.01)

  X_new <- matrix(rnorm(10 * p), 10, p)
  colnames(X_new) <- paste0("V", 1:p)

  preds <- predict(fit, newdata = X_new)

  expect_true(is.matrix(preds))
  expect_equal(dim(preds), c(10, n_tasks))
  expect_equal(colnames(preds), paste0("T", 1:n_tasks))
})

test_that("predict.orthoMTL aligns columns correctly", {
  set.seed(42)
  n <- 20; p <- 4; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- c("A", "B", "C", "D")
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  fit <- orthoMTL(X, Y, lambda = 0.01)

  # Reversed column order
  X_new <- matrix(rnorm(5 * p), 5, p)
  colnames(X_new) <- c("D", "C", "B", "A")

  preds <- predict(fit, newdata = X_new)
  expect_equal(dim(preds), c(5, n_tasks))
})

test_that("predict.orthoMTL errors on missing features", {
  set.seed(42)
  n <- 20; p <- 4; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- c("A", "B", "C", "D")
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  fit <- orthoMTL(X, Y, lambda = 0.01)

  X_new <- matrix(rnorm(5 * 3), 5, 3)
  colnames(X_new) <- c("A", "B", "C")  # missing D

  expect_error(predict(fit, newdata = X_new), "missing")
})

test_that("predict.orthoMTL warns on extra features", {
  set.seed(42)
  n <- 20; p <- 4; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- c("A", "B", "C", "D")
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  fit <- orthoMTL(X, Y, lambda = 0.01)

  X_new <- matrix(rnorm(5 * 5), 5, 5)
  colnames(X_new) <- c("A", "B", "C", "D", "EXTRA")

  expect_warning(predict(fit, newdata = X_new), "not used by the model")
})

test_that("predict.orthoMTL errors on dimension mismatch without names", {
  set.seed(42)
  n <- 20; p <- 4; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)  # no colnames
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)

  fit <- orthoMTL(X, Y, lambda = 0.01)

  X_wrong <- matrix(rnorm(5 * 3), 5, 3)
  expect_error(predict(fit, newdata = X_wrong), "columns")
})

test_that("predict.orthoMTL applies monotonicity in survival mode", {
  set.seed(42)
  n <- 30; p <- 5; n_tasks <- 3
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)

  SurvTime <- rexp(n, 0.1)
  Event <- rbinom(n, 1, 0.7)
  Y <- create_longitudinal_labels(SurvTime, Event, c(4, 6, 10))
  W <- create_indicator_matrix(Y)

  fit <- orthoMTL(X, Y, lambda = 0.01, survival = TRUE, censored.mat = W)

  X_new <- matrix(rnorm(10 * p), 10, p)
  colnames(X_new) <- paste0("V", 1:p)
  preds <- predict(fit, newdata = X_new)

  # Predictions should be non-increasing across columns for each row
  for (i in seq_len(nrow(preds))) {
    diffs <- diff(preds[i, ])
    expect_true(all(diffs <= .Machine$double.eps * 100),
                info = paste("Row", i, "is not non-increasing"))
  }
})
