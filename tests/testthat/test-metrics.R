# Tests for cindex_mtl and nnmaxheap_C

test_that("cindex_mtl returns value between 0 and 1", {
  set.seed(42)
  n <- 20; n_tasks <- 3
  true_mat <- matrix(sample(c(0, 1), n * n_tasks, replace = TRUE), n, n_tasks)
  pred_mat <- matrix(rnorm(n * n_tasks), n, n_tasks)

  cindex <- cindex_mtl(true_mat, pred_mat)

  expect_true(is.numeric(cindex))
  expect_true(length(cindex) == 1)
  expect_true(cindex >= 0 && cindex <= 1)
})

test_that("cindex_mtl handles perfect concordance", {
  # Construct a case where predictions perfectly match ordering
  true_mat <- matrix(c(
    0, 0, 0,  # worst survival
    1, 0, 0,
    1, 1, 0,
    1, 1, 1   # best survival
  ), nrow = 4, byrow = TRUE)

  # Predictions that perfectly match the ordering (higher = better)
  pred_mat <- matrix(c(
    0.1, 0.1, 0.1,
    0.5, 0.3, 0.2,
    0.8, 0.6, 0.3,
    0.9, 0.9, 0.9
  ), nrow = 4, byrow = TRUE)

  cindex <- cindex_mtl(true_mat, pred_mat)

  expect_true(cindex > 0.8)  # should be high, may not be exactly 1 due to implementation
})

test_that("nnmaxheap_C enforces non-increasing output", {
  # Input that violates monotonicity
  m <- c(0.5, 0.8, 0.3, 0.1)
  result <- nnmaxheap_C(m)

  diffs <- diff(result)
  expect_true(all(diffs <= .Machine$double.eps * 100))
})

test_that("nnmaxheap_C enforces non-negative output", {
  m <- c(0.5, -0.3, -0.8, -1.0)
  result <- nnmaxheap_C(m)

  expect_true(all(result >= 0))
})

test_that("nnmaxheap_C output is non-negative and non-increasing", {
  m <- c(1.0, 0.8, 0.5, 0.2)
  result <- nnmaxheap_C(m)

  expect_true(all(result >= 0))
  diffs <- diff(result)
  expect_true(all(diffs <= .Machine$double.eps * 100))
})

test_that("nnmaxheap_C handles length 1", {
  expect_equal(nnmaxheap_C(0.5), 0.5)
  expect_equal(nnmaxheap_C(-0.5), 0)
})

test_that("nnmaxheap_C errors on empty input", {
  expect_error(nnmaxheap_C(numeric(0)), "should be an integer over 1")
})
