# Tests for simulate_mtl

test_that("simulate_mtl returns S3 class with expected structure", {
  sim <- simulate_mtl(n = 50, p = 10, n_signals = 3, seed = 42)

  expect_s3_class(sim, "simulated_mtl")
  expect_true(is.matrix(sim$X))
  expect_true(is.numeric(sim$SurvTime))
  expect_true(is.numeric(sim$Event))
  expect_true(is.numeric(sim$treatment))
  expect_true(is.list(sim$ground_truth))
  expect_equal(sim$n, 50)
  expect_equal(sim$p, 10)
  expect_equal(sim$n_signals, 3)
})

test_that("simulate_mtl X has correct dimensions and names", {
  sim <- simulate_mtl(n = 40, p = 15, n_continuous = 2, seed = 1)

  expect_equal(dim(sim$X), c(40, 16))  # 15 features + treatment
  expect_equal(ncol(sim$X), sim$p + 1)
  expect_equal(colnames(sim$X), sim$feature_names)

  # First 2 columns should be continuous
  expect_true(all(grepl("^cont_", sim$feature_names[1:2])))

  # Last column is treatment
  expect_equal(sim$feature_names[length(sim$feature_names)], "treatment")
})

test_that("simulate_mtl binary features are 0/1", {
  sim <- simulate_mtl(n = 100, p = 20, n_continuous = 1, seed = 42)

  binary_cols <- grep("^mut_", colnames(sim$X))
  for (j in binary_cols) {
    vals <- unique(sim$X[, j])
    expect_true(all(vals %in% c(0, 1)),
                info = paste("Column", colnames(sim$X)[j],
                             "has non-binary values"))
  }
})

test_that("simulate_mtl treatment is balanced", {
  sim <- simulate_mtl(n = 200, p = 10, seed = 42)

  trt_table <- table(sim$treatment)
  expect_equal(as.numeric(trt_table), c(100, 100))
})

test_that("simulate_mtl SurvTime is positive and Event is binary", {
  sim <- simulate_mtl(n = 100, p = 10, seed = 42)

  expect_true(all(sim$SurvTime > 0))
  expect_true(all(sim$Event %in% c(0, 1)))
})

test_that("simulate_mtl censoring rate is reasonable", {
  sim <- simulate_mtl(n = 500, p = 10, censoring_max = 25, seed = 42)

  event_rate <- mean(sim$Event)
  # Should be between 20% and 90% for reasonable settings
  expect_true(event_rate > 0.2 && event_rate < 0.9,
              info = paste("Event rate:", event_rate))
})

test_that("simulate_mtl ground truth has correct structure", {
  sim <- simulate_mtl(n = 50, p = 20, n_signals = 5, seed = 42)

  gt <- sim$ground_truth

  expect_true(is.matrix(gt$coefficients))
  expect_equal(dim(gt$coefficients), c(21, length(sim$thresholds)))
  expect_equal(length(gt$signal_features), 5)
  expect_equal(length(gt$null_features), 15)
  expect_equal(length(gt$effect_types), 5)

  # Signal features should have non-zero coefficients
  for (feat in gt$signal_features) {
    expect_true(any(gt$coefficients[feat, ] != 0),
                info = paste(feat, "has all-zero coefficients"))
  }

  # Null features should have all-zero coefficients
  for (feat in gt$null_features) {
    expect_true(all(gt$coefficients[feat, ] == 0),
                info = paste(feat, "has non-zero coefficients"))
  }

  # Treatment should have non-zero coefficients
  expect_true(any(gt$coefficients["treatment", ] != 0))
})

test_that("simulate_mtl effect types cycle through templates", {
  sim <- simulate_mtl(n = 50, p = 20, n_signals = 10, seed = 42)

  types <- sim$ground_truth$effect_types
  valid_types <- c("early", "late", "constant", "switch", "decreasing")

  expect_true(all(types %in% valid_types))

  # With 10 signals and 5 templates, each template should appear twice
  type_counts <- table(types)
  expect_true(all(type_counts == 2))
})

test_that("simulate_mtl is reproducible with same seed", {
  sim1 <- simulate_mtl(n = 50, p = 10, seed = 123)
  sim2 <- simulate_mtl(n = 50, p = 10, seed = 123)

  expect_identical(sim1$X, sim2$X)
  expect_identical(sim1$SurvTime, sim2$SurvTime)
  expect_identical(sim1$Event, sim2$Event)
  expect_identical(sim1$ground_truth$coefficients,
                   sim2$ground_truth$coefficients)
})

test_that("simulate_mtl produces different data with different seed", {
  sim1 <- simulate_mtl(n = 50, p = 10, seed = 1)
  sim2 <- simulate_mtl(n = 50, p = 10, seed = 2)

  expect_false(identical(sim1$X, sim2$X))
})

test_that("simulate_mtl validates inputs", {
  expect_error(simulate_mtl(n = 5), "at least 10")
  expect_error(simulate_mtl(p = 0), "at least 1")
  expect_error(simulate_mtl(n_signals = 50, p = 10), "cannot exceed")
  expect_error(simulate_mtl(thresholds = c(10, 5)), "strictly increasing")
  expect_error(simulate_mtl(thresholds = c(5)), "At least 2")
  expect_error(simulate_mtl(baseline_hazard = -1), "positive")
  expect_error(simulate_mtl(censoring_max = 0), "positive")
})

test_that("simulate_mtl integrates with orthoMTL workflow", {
  sim <- simulate_mtl(n = 50, p = 10, n_signals = 3, seed = 42)

  Y <- create_longitudinal_labels(sim$SurvTime, sim$Event, sim$thresholds)
  W <- create_indicator_matrix(Y)
  K <- create_constraint_matrix(length(sim$thresholds))

  fit <- orthoMTL(sim$X, Y, lambda = 1e-3, K = K,
                  survival = TRUE, censored.mat = W)

  expect_s3_class(fit, "orthoMTL")
  expect_equal(dim(fit$B), c(ncol(sim$X), length(sim$thresholds)))
})

test_that("print.simulated_mtl produces output", {
  sim <- simulate_mtl(n = 50, p = 10, n_signals = 3, seed = 42)

  expect_output(print(sim), "Simulated multi-task survival data")
  expect_output(print(sim), "Signal features")
  expect_output(print(sim), "Event rate")
  expect_output(print(sim), "Treatment effect")
})
