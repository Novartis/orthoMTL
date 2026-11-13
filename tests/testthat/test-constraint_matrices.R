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

# Tests for validate_K

test_that("validate_K: identity is convex with rho = Inf", {
  chk <- validate_K(diag(4))
  expect_true(chk$convex)
  expect_equal(chk$min_eigenvalue, 1)
  expect_equal(chk$lambda_max_offdiag, 0)
  expect_equal(chk$rho, Inf)
})

test_that("validate_K: uniform K has threshold T - 1 (paper's delta / (T - 1))", {
  T <- 6
  for (d in c(T - 1.5, T - 1, T + 1)) {
    K <- matrix(1, T, T)
    diag(K) <- d
    chk <- validate_K(K)
    expect_equal(chk$lambda_max_offdiag, T - 1)
    expect_equal(chk$rho, d / (T - 1))
    expect_equal(chk$convex, d >= T - 1)
  }
})

test_that("validate_K: convex boundary at diag = lambda_max_offdiag", {
  chk0 <- validate_K(create_constraint_matrix(5))
  expect_false(chk0$convex)
  expect_lt(chk0$rho, 1)

  K <- create_constraint_matrix(5, diag_val = chk0$lambda_max_offdiag)
  chk1 <- validate_K(K)
  expect_true(chk1$convex)
  expect_equal(chk1$rho, 1)
  expect_equal(chk1$min_eigenvalue, 0, tolerance = 1e-8)
})

test_that("validate_K: non-uniform diagonal uses the companion matrix", {
  # Two tasks: convex iff d1 * d2 >= k^2 (companion matrix is PSD)
  K <- matrix(c(4, 2, 2, 1), 2, 2)          # d1 * d2 = 4 = k^2 -> boundary
  expect_true(validate_K(K)$convex)
  K[1, 1] <- 3                               # 3 < 4 -> non-convex
  expect_false(validate_K(K)$convex)
})

test_that("validate_K: rejects invalid input", {
  expect_error(validate_K(1:4), "square numeric matrix")
  expect_error(validate_K(matrix(1, 2, 3)), "square numeric matrix")
  expect_error(validate_K(matrix(c(1, NA, NA, 1), 2)), "NA")
  expect_error(validate_K(matrix(c(1, 2, 0, 1), 2)), "symmetric")
  expect_error(validate_K(matrix(c(1, -1, -1, 1), 2)), "non-negative")
})
