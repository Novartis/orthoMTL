# Tests for plotting functions

test_that("plot_heatmap returns ggplot from orthoMTL object", {
  set.seed(42)
  X <- matrix(rnorm(60), 12, 5)
  colnames(X) <- paste0("V", 1:5)
  Y <- matrix(rnorm(36), 12, 3)
  colnames(Y) <- paste0("T", 1:3)
  fit <- orthoMTL(X, Y, lambda = 0.01)

  p <- plot_heatmap(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap returns ggplot from raw matrix", {
  Bt <- matrix(rnorm(15), 5, 3)
  rownames(Bt) <- paste0("V", 1:5)
  colnames(Bt) <- paste0("T", 1:3)

  p <- plot_heatmap(Bt)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap respects reorder = FALSE", {
  Bt <- matrix(c(3, 1, 2, 6, 4, 5), 3, 2)
  rownames(Bt) <- c("A", "B", "C")
  colnames(Bt) <- c("T1", "T2")

  p_ordered <- plot_heatmap(Bt, reorder = TRUE)
  p_raw <- plot_heatmap(Bt, reorder = FALSE)

  # Both should be valid ggplots
  expect_s3_class(p_ordered, "ggplot")
  expect_s3_class(p_raw, "ggplot")

  # Extract row order from the data
  raw_levels <- levels(p_raw$data$X)
  expect_equal(raw_levels, c("A", "B", "C"))
})

test_that("plot_heatmap generates names when absent", {
  Bt <- matrix(rnorm(12), 4, 3)  # no row/col names

  p <- plot_heatmap(Bt)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap errors on invalid input", {
  expect_error(plot_heatmap("not_a_matrix"), "must be")
})

test_that("plot_correlation returns ggplot", {
  set.seed(42)
  X <- matrix(rnorm(60), 12, 5)
  colnames(X) <- paste0("V", 1:5)
  Y <- matrix(rnorm(48), 12, 4)
  colnames(Y) <- paste0("T", 1:4)
  fit <- orthoMTL(X, Y, lambda = 0.01)

  p <- plot_correlation(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_correlation accepts custom midpoint and limits", {
  Bt <- matrix(rnorm(20), 5, 4)
  colnames(Bt) <- paste0("T", 1:4)

  p <- plot_correlation(Bt, midpoint = 0.3, limits = c(-0.5, 1.5))
  expect_s3_class(p, "ggplot")
})

test_that("plot_correlation errors on invalid input", {
  expect_error(plot_correlation(list(a = 1)), "must be")
})

test_that("plot_prediction returns ggplot", {
  Mt <- matrix(rnorm(30), 10, 3)
  rownames(Mt) <- paste0("pt_", 1:10)
  colnames(Mt) <- c("T1", "T2", "T3")

  p <- plot_prediction(Mt)
  expect_s3_class(p, "ggplot")
})

test_that("plot_prediction generates names when absent", {
  Mt <- matrix(rnorm(12), 4, 3)  # no names

  p <- plot_prediction(Mt)
  expect_s3_class(p, "ggplot")
})

test_that("plot_prediction errors on non-matrix", {
  expect_error(plot_prediction(data.frame(a = 1)), "must be a numeric matrix")
})

test_that("plot_bootstrap returns list of ggplots from S3 object", {
  set.seed(42)
  n <- 20; p <- 3; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  colnames(Y) <- paste0("T", 1:n_tasks)

  boot_res <- bootstrap_orthoMTL(
    X = X, Y = Y, lambda = 0.01,
    n_repeats = 3, n_cores = 1, verbose = FALSE
  )

  plots <- plot_bootstrap(boot_res)
  expect_true(is.list(plots))
  expect_true(length(plots) >= 1)
  expect_s3_class(plots[[1]], "ggplot")
})

test_that("plot_bootstrap accepts raw data.frame", {
  df <- data.frame(
    id = rep(c("V1", "V2"), each = 8),
    time = rep(rep(c("T1", "T2"), each = 2), 4),
    coeff = rnorm(16),
    group = rep(c("real", "null"), 8),
    stringsAsFactors = FALSE
  )

  plots <- plot_bootstrap(df)
  expect_true(is.list(plots))
  expect_s3_class(plots[[1]], "ggplot")
})

test_that("plot_bootstrap subsets features correctly", {
  set.seed(42)
  n <- 20; p <- 5; n_tasks <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  Y <- matrix(rnorm(n * n_tasks), n, n_tasks)
  colnames(Y) <- paste0("T", 1:n_tasks)

  boot_res <- bootstrap_orthoMTL(
    X = X, Y = Y, lambda = 0.01,
    n_repeats = 3, n_cores = 1, verbose = FALSE
  )

  plots <- plot_bootstrap(boot_res, features = c("V1", "V3"))
  expect_true(is.list(plots))
  expect_s3_class(plots[[1]], "ggplot")
})

test_that("plot_bootstrap warns on missing features", {
  df <- data.frame(
    id = rep("V1", 4),
    time = rep(c("T1", "T2"), 2),
    coeff = rnorm(4),
    group = rep(c("real", "null"), 2),
    stringsAsFactors = FALSE
  )

  expect_warning(
    plot_bootstrap(df, features = c("V1", "NONEXISTENT")),
    "not found"
  )
})

test_that("plot_bootstrap respects batch_size", {
  df <- data.frame(
    id = rep(paste0("V", 1:10), each = 4),
    time = rep(rep(c("T1", "T2"), each = 2), 10),
    coeff = rnorm(40),
    group = rep(c("real", "null"), 20),
    stringsAsFactors = FALSE
  )

  plots <- suppressWarnings(plot_bootstrap(df, batch_size = 4))
  expect_equal(length(plots), 3)
})

test_that("plot_bootstrap errors on invalid input", {
  expect_error(plot_bootstrap("not_valid"), "must be")
})

test_that("plot_bootstrap errors on data.frame with missing columns", {
  df <- data.frame(id = "V1", time = "T1", coeff = 1.0)
  # missing 'group' column
  expect_error(plot_bootstrap(df), "missing required columns")
})
