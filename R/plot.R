#' Coefficient Heatmap for orthoMTL
#'
#' Displays the coefficient matrix as a heatmap. Features (rows) can
#' optionally be reordered by their mean coefficient across tasks.
#'
#' @param x Either a fitted \code{"orthoMTL"} object or a numeric
#'   coefficient matrix with dimensions \code{p x numTasks}. If an
#'   \code{"orthoMTL"} object, coefficients are extracted via
#'   \code{coef()}.
#' @param reorder Logical. If \code{TRUE} (default), rows are sorted
#'   by mean coefficient across tasks (ascending).
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{coef.orthoMTL}}, \code{\link{plot_correlation}}
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 100; p <- 10; n_tasks <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' Y <- X %*% matrix(rnorm(p * n_tasks), p) + matrix(rnorm(n * n_tasks), n) * 0.1
#' colnames(Y) <- c("T1", "T2", "T3")
#' K <- matrix(1, n_tasks, n_tasks); diag(K) <- 0.5
#' fit <- orthoMTL(X, Y, lambda = 1e-3, K = K)
#'
#' # From fitted object
#' plot_heatmap(fit)
#'
#' # From raw matrix
#' plot_heatmap(coef(fit), reorder = FALSE)
plot_heatmap <- function(x, reorder = TRUE) {

  # Extract coefficient matrix
  if (inherits(x, "orthoMTL")) {
    Bt <- coef(x)
  } else if (is.matrix(x)) {
    Bt <- x
  } else {
    stop("'x' must be an 'orthoMTL' object or a numeric matrix.",
         call. = FALSE)
  }

  # Ensure names exist
  if (is.null(rownames(Bt))) {
    rownames(Bt) <- paste0("V", seq_len(nrow(Bt)))
  }
  if (is.null(colnames(Bt))) {
    colnames(Bt) <- paste0("T", seq_len(ncol(Bt)))
  }

  # Optionally reorder rows by mean coefficient
  if (reorder) {
    row_order <- order(apply(Bt, 1, mean), decreasing = FALSE)
    Bt <- Bt[row_order, , drop = FALSE]
  }

  # Build long-format data
  data <- expand.grid(
    X = factor(rownames(Bt), levels = rownames(Bt)),
    Y = factor(colnames(Bt), levels = colnames(Bt))
  )
  data$Z <- as.vector(Bt)

  ggplot2::ggplot(data, ggplot2::aes(.data$Y, .data$X, fill = .data$Z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "red", mid = "white", high = "blue"
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = "weight") +
    ggplot2::theme_minimal()
}


#' Task Correlation Map for orthoMTL
#'
#' Displays pairwise distances between task coefficient vectors as a
#' heatmap. Distance is computed as \code{1 - cosine_similarity}.
#' Values near 0 indicate similar coefficient profiles; values near 1
#' indicate orthogonal profiles.
#'
#' @param x Either a fitted \code{"orthoMTL"} object or a numeric
#'   coefficient matrix with dimensions \code{p x numTasks}.
#' @param midpoint Midpoint for the colour scale. Default: \code{0.5}.
#' @param limits Numeric vector of length 2 for the colour scale limits.
#'   Default: \code{c(0, 1)}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{plot_heatmap}}
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 100; p <- 10; n_tasks <- 4
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' Y <- X %*% matrix(rnorm(p * n_tasks), p) + matrix(rnorm(n * n_tasks), n) * 0.1
#' colnames(Y) <- c("T1", "T2", "T3", "T4")
#' K <- matrix(1, n_tasks, n_tasks); diag(K) <- 0.5
#' fit <- orthoMTL(X, Y, lambda = 1e-3, K = K)
#'
#' plot_correlation(fit)
plot_correlation <- function(x, midpoint = 0.5, limits = c(0, 1)) {

  # Extract coefficient matrix
  if (inherits(x, "orthoMTL")) {
    Bt <- coef(x)
  } else if (is.matrix(x)) {
    Bt <- x
  } else {
    stop("'x' must be an 'orthoMTL' object or a numeric matrix.",
         call. = FALSE)
  }

  if (is.null(colnames(Bt))) {
    colnames(Bt) <- paste0("T", seq_len(ncol(Bt)))
  }

  # Normalise columns to unit length
  norms <- apply(Bt, 2, function(col) sqrt(sum(col^2)))
  # Guard against zero-norm columns
  norms[norms == 0] <- 1
  normBt <- t(t(Bt) / norms)

  # Cosine distance: 1 - cosine_similarity
  cos_sim <- crossprod(normBt)
  dist_mat <- 1 - cos_sim

  # Clamp to [0, 1] to prevent floating-point overflow → NA in ggplot
  dist_mat[dist_mat < 0] <- 0
  dist_mat[dist_mat > 1] <- 1

  # Reverse row order for display
  task_names <- colnames(Bt)
  task_names_rev <- rev(task_names)
  dist_mat_rev <- dist_mat[task_names_rev, , drop = FALSE]

  # Build long-format data
  data <- expand.grid(
    X = factor(task_names_rev, levels = task_names_rev),
    Y = factor(task_names, levels = task_names)
  )
  data$distance <- as.vector(dist_mat_rev)

  ggplot2::ggplot(data, ggplot2::aes(.data$Y, .data$X, fill = .data$distance)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "red", mid = "white", high = "blue",
      midpoint = midpoint,
      limits = limits
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = "distance") +
    ggplot2::theme_minimal()
}


#' Prediction Swimmer Plot for orthoMTL
#'
#' Displays a prediction matrix as a heatmap (swimmer plot), where
#' rows are patients and columns are tasks/thresholds.
#'
#' @param x A numeric prediction matrix with dimensions
#'   \code{n x numTasks}, typically produced by
#'   \code{\link{predict.orthoMTL}}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{predict.orthoMTL}}
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 50; p <- 10; n_tasks <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' Y <- X %*% matrix(rnorm(p * n_tasks), p) + matrix(rnorm(n * n_tasks), n) * 0.1
#' colnames(Y) <- c("T1", "T2", "T3")
#' K <- matrix(1, n_tasks, n_tasks); diag(K) <- 0.5
#' fit <- orthoMTL(X, Y, lambda = 1e-3, K = K)
#'
#' preds <- predict(fit, newdata = X)
#' plot_prediction(preds)
plot_prediction <- function(x) {

  if (!is.matrix(x)) {
    stop("'x' must be a numeric matrix.", call. = FALSE)
  }

  # Ensure names exist
  if (is.null(rownames(x))) {
    rownames(x) <- paste0("obs_", seq_len(nrow(x)))
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("T", seq_len(ncol(x)))
  }

  # Build long-format data
  data <- expand.grid(
    X = factor(rownames(x), levels = rownames(x)),
    Y = factor(colnames(x), levels = colnames(x))
  )
  data$Z <- as.vector(x)

  ggplot2::ggplot(data, ggplot2::aes(.data$Y, .data$X, fill = .data$Z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "red", mid = "white", high = "blue"
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = "predicted") +
    ggplot2::theme_minimal()
}


#' Bootstrap Coefficient Comparison Plot
#'
#' Displays faceted line plots comparing real (bootstrapped) vs null
#' (permuted) coefficient trajectories across tasks for selected or
#' all features. Each panel shows the mean and standard error of the
#' coefficient at each task/threshold.
#'
#' @param x Either a \code{"bootstrap_orthoMTL"} object (as returned
#'   by \code{\link{bootstrap_orthoMTL}}) or a \code{data.frame} with
#'   columns \code{id}, \code{time}, \code{coeff}, and \code{group}.
#' @param features A character vector of feature names to plot. If
#'   \code{NULL} (default), all features are plotted.
#' @param batch_size Number of features per facet page. Default:
#'   \code{9}.
#'
#' @return A list of \code{ggplot2} objects, one per batch/page.
#'
#' @details
#' Blue lines show coefficients from models fitted on bootstrapped
#' (resampled) data — reflecting estimation variability under real
#' signal. Grey lines show coefficients from models fitted on
#' permuted outcomes — reflecting the null distribution.
#'
#' When the blue (real) band is clearly separated from the grey
#' (null) band, the feature's coefficient is distinguishable from
#' noise at that task/threshold.
#'
#' @seealso \code{\link{bootstrap_orthoMTL}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 30; p <- 5; n_tasks <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' SurvTime <- rexp(n, rate = 0.1)
#' Event <- rbinom(n, 1, 0.7)
#' thresholds <- c(4, 6, 10)
#' Y <- create_longitudinal_labels(SurvTime, Event, thresholds)
#' W <- create_indicator_matrix(Y)
#' K <- create_constraint_matrix(n_tasks)
#' boot_res <- bootstrap_orthoMTL(
#'   X = X, Y = Y, lambda = 1e-3, step_size = 0.1,
#'   K = K, survival = TRUE, censored.mat = W,
#'   n_repeats = 5, n_cores = 1, verbose = FALSE
#' )
#' plots <- plot_bootstrap(boot_res)
#' plots[[1]]
#' }
plot_bootstrap <- function(x, features = NULL, batch_size = 9) {

  # Extract data.frame from S3 object or use directly
  if (inherits(x, "bootstrap_orthoMTL")) {
    df <- x$results
  } else if (is.data.frame(x)) {
    required_cols <- c("id", "time", "coeff", "group")
    missing_cols <- setdiff(required_cols, colnames(x))
    if (length(missing_cols) > 0) {
      stop("data.frame is missing required columns: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    df <- x
  } else {
    stop("'x' must be a 'bootstrap_orthoMTL' object or a data.frame ",
         "with columns: id, time, coeff, group.", call. = FALSE)
  }

  # Subset features if requested
  if (!is.null(features)) {
    missing_feats <- setdiff(features, unique(df$id))
    if (length(missing_feats) > 0) {
      warning("The following features were not found and will be skipped: ",
              paste(missing_feats, collapse = ", "), call. = FALSE)
    }
    df <- df[df$id %in% features, , drop = FALSE]
    feature_order <- features[features %in% unique(df$id)]
  } else {
    feature_order <- unique(df$id)
  }

  n_features <- length(feature_order)

  if (n_features == 0) {
    warning("No features to plot.", call. = FALSE)
    return(list())
  }

  # Batch features into pages
  n_batches <- ceiling(n_features / batch_size)
  plot_list <- vector("list", n_batches)

  for (b in seq_len(n_batches)) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, n_features)
    batch_features <- feature_order[start_idx:end_idx]

    batch_df <- df[df$id %in% batch_features, , drop = FALSE]
    batch_df$feature <- factor(batch_df$id, levels = batch_features)
    # Safely convert: only reorder numerically if all values can be parsed as numbers
    times_as_num <- suppressWarnings(as.numeric(as.character(batch_df[["time"]])))
    if (!any(is.na(times_as_num))) {
      # All numeric — sort numerically
      batch_df[["time"]] <- factor(batch_df[["time"]],
                                   levels = sort(unique(as.character(batch_df[["time"]])[order(times_as_num)])))
    } else {
      # Non-numeric labels — preserve original order
      batch_df[["time"]] <- factor(batch_df[["time"]],
                                   levels = unique(batch_df[["time"]]))
    }

    plot_list[[b]] <- ggplot2::ggplot(
      batch_df,
      ggplot2::aes(
        x = .data$time,
        y = .data$coeff,
        colour = .data$group,
        group = .data$group
      )
    ) +
      ggplot2::stat_summary(fun = mean, geom = "line") +
      ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                            width = 0.2) +
      ggplot2::facet_wrap(~ feature, scales = "free_y") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::xlab("") +
      ggplot2::ylab("model coefficient") +
      ggplot2::theme_minimal()
  }

  return(plot_list)
}
