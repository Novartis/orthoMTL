# Tests for create_longitudinal_labels and create_indicator_matrix

test_that("create_longitudinal_labels produces correct dimensions", {
  SurvTime <- c(2, 5, 8, 12, 20)
  Event <- c(1, 1, 0, 1, 0)
  thresholds <- c(4, 6, 10)

  Y <- create_longitudinal_labels(SurvTime, Event, thresholds)

  expect_true(is.matrix(Y))
  expect_equal(dim(Y), c(5, 3))
  expect_equal(colnames(Y), as.character(thresholds))
})

test_that("create_longitudinal_labels encodes correctly", {
  SurvTime <- c(2, 5, 8, 12)
  Event <-    c(1, 1, 0, 1)
  thresholds <- c(4, 6, 10)

  Y <- create_longitudinal_labels(SurvTime, Event, thresholds)

  # Patient 1: event at 2, before all thresholds

  expect_equal(unname(Y[1, ]), c(0, 0, 0))

  # Patient 2: event at 5, before threshold 6 and 10
  expect_equal(unname(Y[2, ]), c(1, 0, 0))

  # Patient 3: censored at 8, before threshold 10
  expect_equal(unname(Y[3, ]), c(1, 1, NA))

  # Patient 4: event at 12, after all thresholds
  expect_equal(unname(Y[4, ]), c(1, 1, 1))
})

test_that("create_longitudinal_labels validates inputs", {
  expect_error(
    create_longitudinal_labels(c(-1, 5), c(1, 0), c(4, 6)),
    "positive"
  )
  expect_error(
    create_longitudinal_labels(c(2, 5), c(1, 2), c(4, 6)),
    "binary"
  )
})

test_that("create_indicator_matrix marks NAs as zero", {
  Y <- matrix(c(1, 0, NA, 1, NA, 1), nrow = 3, ncol = 2)
  W <- create_indicator_matrix(Y)

  expect_equal(dim(W), dim(Y))
  expect_equal(W[1, ], c(1, 1))
  expect_equal(W[2, ], c(1, 0))
  expect_equal(W[3, ], c(0, 1))
})

test_that("create_indicator_matrix handles all-observed data", {
  Y <- matrix(c(1, 0, 1, 0), nrow = 2, ncol = 2)
  W <- create_indicator_matrix(Y)

  expect_true(all(W == 1))
})
