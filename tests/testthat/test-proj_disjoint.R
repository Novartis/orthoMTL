# Tests for proj_disjoint (internal function)

test_that("proj_disjoint returns list with w and v", {
  w <- matrix(rnorm(12), 4, 3)
  v <- matrix(abs(rnorm(12)), 4, 3)

  result <- orthoMTL:::proj_disjoint(w, v)

  expect_true(is.list(result))
  expect_true("w" %in% names(result))
  expect_true("v" %in% names(result))
  expect_equal(dim(result$w), dim(w))
  expect_equal(dim(result$v), dim(v))
})

test_that("proj_disjoint projects negative v to zero", {
  w <- matrix(c(1, 2, 3, 4), 2, 2)
  v <- matrix(c(-1, 2, -3, 4), 2, 2)

  result <- orthoMTL:::proj_disjoint(w, v)

  # Where v was negative, both should be 0
  expect_equal(result$v[1, 1], 0)
  expect_equal(result$w[1, 1], 0)
  expect_equal(result$v[1, 2], 0)
  expect_equal(result$w[1, 2], 0)
})

test_that("proj_disjoint v is non-negative after projection", {
  set.seed(42)
  w <- matrix(rnorm(20), 5, 4)
  v <- matrix(rnorm(20), 5, 4)

  result <- orthoMTL:::proj_disjoint(w, v)

  expect_true(all(result$v >= 0))
})
