# Tests for simulate_mtl mode = regression / classification (SIM-01)

test_that("default mode is survival and back-compatible", {
  sim <- simulate_mtl(n = 50, p = 10, n_signals = 3, seed = 42)
  expect_equal(sim$mode, "survival")
  expect_true(is.numeric(sim$SurvTime))
  expect_true(is.numeric(sim$Event))
  expect_null(sim$Y)
})

test_that("regression mode returns a numeric Y and no survival fields", {
  sim <- simulate_mtl(n = 60, p = 10, n_signals = 3,
                      mode = "regression", seed = 42)
  expect_equal(sim$mode, "regression")
  expect_true(is.matrix(sim$Y))
  expect_equal(dim(sim$Y), c(60, length(sim$thresholds)))
  expect_null(sim$SurvTime)
  expect_null(sim$Event)
  expect_true(all(is.finite(sim$Y)))
})

test_that("regression noise_sd influences residual spread", {
  s_lo <- simulate_mtl(n = 200, p = 8, mode = "regression",
                       noise_sd = 0.1, seed = 1)
  s_hi <- simulate_mtl(n = 200, p = 8, mode = "regression",
                       noise_sd = 5, seed = 1)
  signal_lo <- s_lo$X %*% s_lo$ground_truth$coefficients
  signal_hi <- s_hi$X %*% s_hi$ground_truth$coefficients
  expect_true(sd(s_hi$Y - signal_hi) > sd(s_lo$Y - signal_lo))
})

test_that("classification mode returns labels in {-1, +1}", {
  sim <- simulate_mtl(n = 80, p = 10, n_signals = 3,
                      mode = "classification", seed = 42)
  expect_equal(sim$mode, "classification")
  expect_true(is.matrix(sim$Y))
  expect_equal(dim(sim$Y), c(80, length(sim$thresholds)))
  expect_true(all(sim$Y %in% c(-1, 1)))
  expect_null(sim$SurvTime)
})

test_that("non-survival modes keep the same ground-truth machinery", {
  sim <- simulate_mtl(n = 50, p = 20, n_signals = 5,
                      mode = "regression", seed = 42)
  gt <- sim$ground_truth
  expect_equal(dim(gt$coefficients), c(21, length(sim$thresholds)))
  expect_equal(length(gt$signal_features), 5)
})

test_that("regression mode round-trips through orthoMTL + metrics", {
  sim <- simulate_mtl(n = 80, p = 10, n_signals = 4,
                      mode = "regression", noise_sd = 0.2, seed = 42)
  fit <- orthoMTL(sim$X, sim$Y, lambda = 1e-3)
  pred <- predict(fit, newdata = sim$X)
  expect_true(rmse_mtl(sim$Y, pred) < sd(sim$Y))  # better than trivial
})

test_that("classification mode round-trips through orthoMTL + auc", {
  sim <- simulate_mtl(n = 120, p = 8, n_signals = 4,
                      mode = "classification", seed = 7)
  fit <- orthoMTL(sim$X, sim$Y, lambda = 1e-3, logistic = TRUE)
  pred <- predict(fit, newdata = sim$X, type = "link")
  expect_true(auc_mtl(sim$Y, pred) > 0.5)  # learns something
})

test_that("print.simulated_mtl adapts to mode", {
  sim_r <- simulate_mtl(n = 40, p = 8, mode = "regression", seed = 1)
  expect_output(print(sim_r), "Simulated multi-task regression data")
  expect_output(print(sim_r), "Response range")

  sim_c <- simulate_mtl(n = 40, p = 8, mode = "classification", seed = 1)
  expect_output(print(sim_c), "Simulated multi-task classification data")
  expect_output(print(sim_c), "Positive-label rate")
})
