# Tests for general-purpose metrics (GEN-01):
# rmse_mtl, r2_mtl, accuracy_mtl, auc_mtl

test_that("rmse_mtl is zero for perfect predictions and positive otherwise", {
  set.seed(1)
  Y <- matrix(rnorm(30), 10, 3)
  expect_equal(rmse_mtl(Y, Y), 0)
  expect_true(rmse_mtl(Y, Y + 0.5) > 0)
})

test_that("rmse_mtl matches the manual pooled formula", {
  Y <- matrix(c(1, 2, 3, 4), 2, 2)
  P <- matrix(c(1, 2, 3, 5), 2, 2)   # one cell off by 1
  expect_equal(rmse_mtl(Y, P), sqrt(1 / 4))
})

test_that("r2_mtl is 1 for perfect fit and <= 1 otherwise", {
  set.seed(1)
  Y <- matrix(rnorm(30), 10, 3)
  expect_equal(r2_mtl(Y, Y), 1)
  P <- Y + matrix(rnorm(30, sd = 0.1), 10, 3)
  r2 <- r2_mtl(Y, P)
  expect_true(r2 < 1 && r2 > 0.9)
})

test_that("r2_mtl errors when true values have zero variance", {
  Y <- matrix(3, 5, 2)
  expect_error(r2_mtl(Y, Y + 1), "zero variance")
})

test_that("accuracy_mtl is 1 when signs all match", {
  Y <- matrix(c(-1, 1, 1, -1), 2, 2)
  P <- Y * 0.3   # same signs, smaller magnitude
  expect_equal(accuracy_mtl(Y, P), 1)
})

test_that("accuracy_mtl works for 0/1 labels and a 0.5 threshold", {
  Y <- matrix(c(0, 1, 1, 0), 2, 2)          # label > 0 == positive
  probs <- matrix(c(0.2, 0.8, 0.9, 0.1), 2, 2)
  expect_equal(accuracy_mtl(Y, probs, threshold = 0.5), 1)
})

test_that("auc_mtl is 1 for perfectly separating scores and 0.5 at chance", {
  Y <- matrix(c(-1, -1, 1, 1), 2, 2)
  P <- matrix(c(0.1, 0.2, 0.8, 0.9), 2, 2)  # positives score higher
  expect_equal(auc_mtl(Y, P), 1)
})

test_that("auc_mtl is invariant to monotone transforms (link vs sigmoid)", {
  set.seed(7)
  Y <- matrix(sample(c(-1, 1), 40, replace = TRUE), 10, 4)
  link <- Y + matrix(rnorm(40), 10, 4)
  resp <- 1 / (1 + exp(-link))
  expect_equal(auc_mtl(Y, link), auc_mtl(Y, resp))
})

test_that("auc_mtl errors when only one class is present", {
  Y <- matrix(1, 4, 2)
  expect_error(auc_mtl(Y, matrix(rnorm(8), 4, 2)), "positive .* and negative")
})

test_that("general metrics ignore NA cells", {
  Y <- matrix(c(1, NA, 3, 4), 2, 2)
  P <- matrix(c(1, 99, 3, 4), 2, 2)   # the NA-aligned cell is wrong but masked
  expect_equal(rmse_mtl(Y, P), 0)
})

test_that("metrics error on dimension mismatch", {
  expect_error(rmse_mtl(matrix(1, 2, 2), matrix(1, 2, 3)),
               "identical dimensions")
})
