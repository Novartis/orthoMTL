# Tests for the orthoMTL `schedule` argument (S-02)

make_reg <- function(seed = 5, n = 120, p = 10, n_tasks = 3) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  B <- matrix(rnorm(p * n_tasks), p, n_tasks)
  Y <- X %*% B + matrix(rnorm(n * n_tasks), n, n_tasks) * 0.3
  list(X = X, Y = Y)
}

test_that("default schedule is 'sqrt' and recorded in hyperparameters", {
  d <- make_reg()
  fit <- orthoMTL(d$X, d$Y, lambda = 1e-2)
  expect_equal(fit$hyperparameters$schedule, "sqrt")
})

test_that("explicit schedule = 'sqrt' equals the default (back-compat)", {
  d <- make_reg()
  fit_default <- orthoMTL(d$X, d$Y, lambda = 1e-2, seed = 1)
  fit_sqrt    <- orthoMTL(d$X, d$Y, lambda = 1e-2, seed = 1, schedule = "sqrt")
  expect_identical(fit_default$B, fit_sqrt$B)
})

test_that("all schedules produce finite coefficients of the right shape", {
  d <- make_reg()
  for (sc in c("sqrt", "log", "const", "linear")) {
    # 'linear' (1/i) may not converge within the budget; that emits the
    # documented max_iter warning, which is expected here -- we only check
    # the returned coefficients are well-formed.
    fit <- suppressWarnings(
      orthoMTL(d$X, d$Y, lambda = 1e-2, schedule = sc,
               stop_no_improve = 200, max_iter = 5e4))
    expect_true(all(is.finite(fit$B)), info = sc)
    expect_equal(dim(fit$B), c(ncol(d$X), ncol(d$Y)), info = sc)
  }
})

test_that("'log' and 'const' reach the sqrt optimum in fewer iterations", {
  d <- make_reg()
  fit_sqrt  <- orthoMTL(d$X, d$Y, lambda = 1e-2, schedule = "sqrt",
                        stop_no_improve = 300, max_iter = 1e5)
  fit_log   <- orthoMTL(d$X, d$Y, lambda = 1e-2, schedule = "log",
                        stop_no_improve = 300, max_iter = 1e5)
  fit_const <- orthoMTL(d$X, d$Y, lambda = 1e-2, schedule = "const",
                        stop_no_improve = 300, max_iter = 1e5)

  # Same optimum (within a small tolerance) ...
  expect_lt(abs(fit_log$obj   - fit_sqrt$obj), 1e-3)
  expect_lt(abs(fit_const$obj - fit_sqrt$obj), 1e-3)
  # ... reached strictly faster than sqrt.
  expect_lt(fit_log$imax,   fit_sqrt$imax)
  expect_lt(fit_const$imax, fit_sqrt$imax)
})

test_that("invalid schedule is rejected", {
  d <- make_reg()
  expect_error(orthoMTL(d$X, d$Y, lambda = 1e-2, schedule = "bogus"),
               "'arg'")
})

test_that("cv_orthoMTL accepts a schedule and rejects a bad one", {
  d <- make_reg(n = 30, p = 4, n_tasks = 2)
  folds <- rep(1:2, length.out = nrow(d$X))
  cv_res <- cv_orthoMTL(
    X.train = d$X, Y.train = d$Y, schedule = "log",
    lambdas = c(0.01), stepsizes = c(0.1), diag_vals = c(1),
    survival = FALSE, folds = folds, n_cores = 1, verbose = FALSE
  )
  expect_s3_class(cv_res, "cv_orthoMTL")
  expect_true(is.finite(cv_res$best$cv_score))

  expect_error(
    cv_orthoMTL(X.train = d$X, Y.train = d$Y, schedule = "bogus",
                lambdas = 0.01, stepsizes = 0.1, diag_vals = 1,
                survival = FALSE, folds = folds, n_cores = 1, verbose = FALSE),
    "'arg'"
  )
})

test_that("bootstrap_orthoMTL accepts a schedule and rejects a bad one", {
  d <- make_reg(n = 30, p = 4, n_tasks = 2)
  colnames(d$X) <- paste0("V", seq_len(ncol(d$X)))
  boot <- bootstrap_orthoMTL(
    X = d$X, Y = d$Y, lambda = 0.01, schedule = "const",
    survival = FALSE, n_repeats = 3, n_cores = 1, verbose = FALSE
  )
  expect_s3_class(boot, "bootstrap_orthoMTL")
  expect_length(boot$obj_real, 3)

  expect_error(
    bootstrap_orthoMTL(X = d$X, Y = d$Y, lambda = 0.01, schedule = "bogus",
                       survival = FALSE, n_repeats = 2, n_cores = 1,
                       verbose = FALSE),
    "'arg'"
  )
})
