# Tests for create_constraint_matrix

test_that("create_constraint_matrix has correct dimensions", {
  K <- create_constraint_matrix(5)

  expect_true(is.matrix(K))
  expect_equal(dim(K), c(5, 5))
})

test_that("create_constraint_matrix is symmetric", {
  K <- create_constraint_matrix(7)

  expect_equal(K, t(K))
})

test_that("create_constraint_matrix diagonal uses diag_val", {
  K1 <- create_constraint_matrix(4, diag_val = 0.5)
  K2 <- create_constraint_matrix(4, diag_val = 3.0)

  expect_equal(diag(K1), rep(0.5, 4))
  expect_equal(diag(K2), rep(3.0, 4))
})

test_that("create_constraint_matrix off-diagonals increase with distance", {
  K <- create_constraint_matrix(5)

  # Adjacent tasks should have smaller off-diagonal than distant tasks
  expect_true(K[1, 2] < K[1, 4])
  expect_true(K[1, 2] < K[1, 5])
  expect_true(K[2, 3] < K[2, 5])
})

test_that("create_constraint_matrix handles 2 tasks", {
  K <- create_constraint_matrix(2)

  expect_equal(dim(K), c(2, 2))
  expect_equal(diag(K), c(0.5, 0.5))
})
