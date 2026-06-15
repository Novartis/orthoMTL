# Tests for predict.orthoMTL type= argument (GEN-03)

make_logistic_fit <- function(seed = 42) {
  set.seed(seed)
  n <- 60; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))
  # labels in {-1, +1}, the encoding the logistic solver expects
  lin <- X[, 1] - X[, 2]
  Y <- matrix(0, n, n_tasks)
  for (k in seq_len(n_tasks)) {
    Y[, k] <- ifelse(stats::runif(n) < 1 / (1 + exp(-lin)), 1, -1)
  }
  list(fit = orthoMTL(X, Y, lambda = 0.01, logistic = TRUE),
       X = X, n_tasks = n_tasks)
}

test_that("type = 'link' is the default and unchanged", {
  obj <- make_logistic_fit()
  p_default <- predict(obj$fit, newdata = obj$X)
  p_link    <- predict(obj$fit, newdata = obj$X, type = "link")
  expect_equal(p_default, p_link)
})

test_that("type = 'response' returns probabilities in (0, 1)", {
  obj <- make_logistic_fit()
  resp <- predict(obj$fit, newdata = obj$X, type = "response")
  expect_true(all(resp > 0 & resp < 1))
  # response is the sigmoid of the link
  link <- predict(obj$fit, newdata = obj$X, type = "link")
  expect_equal(resp, 1 / (1 + exp(-link)))
})

test_that("type = 'class' returns labels in {-1, +1}", {
  obj <- make_logistic_fit()
  cls <- predict(obj$fit, newdata = obj$X, type = "class")
  expect_true(all(cls %in% c(-1, 1)))
  link <- predict(obj$fit, newdata = obj$X, type = "link")
  expect_equal(cls, ifelse(link >= 0, 1, -1))
})

test_that("type = 'class' errors for non-logistic fits", {
  set.seed(1)
  n <- 30; p <- 4; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  fit <- orthoMTL(X, Y, lambda = 0.01)  # regression
  expect_error(predict(fit, newdata = X, type = "class"),
               "only available for logistic")
})

test_that("type = 'response' equals link for non-logistic fits", {
  set.seed(1)
  n <- 30; p <- 4; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  fit <- orthoMTL(X, Y, lambda = 0.01)
  expect_equal(predict(fit, newdata = X, type = "response"),
               predict(fit, newdata = X, type = "link"))
})
