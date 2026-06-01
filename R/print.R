#' Print an orthoMTL Object
#'
#' Displays a compact summary of a fitted \code{orthoMTL} model.
#'
#' @param x A fitted model object of class \code{"orthoMTL"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print orthoMTL
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 100; p <- 10; n_tasks <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X %*% matrix(rnorm(p * n_tasks), p) + matrix(rnorm(n * n_tasks), n) * 0.1
#' fit <- orthoMTL(X, Y, lambda = 1e-3)
#' print(fit)
print.orthoMTL <- function(x, ...) {

  # Determine mode
  mode <- if (isTRUE(x$hyperparameters$survival)) {
    "survival"
  } else if (isTRUE(x$hyperparameters$logistic)) {
    "classification"
  } else {
    "regression"
  }

  cat("orthoMTL model:", x$n_features, "features,",
      x$n_tasks, "tasks (", mode, ")\n", sep = " ")

  conv_str <- if (x$converged) "yes" else "NO"
  cat("Converged:", conv_str, "(", x$imax, "iterations, objective:",
      format(x$obj, digits = 4, nsmall = 4), ")\n", sep = " ")

  pen_str <- paste0("lambda = ", x$hyperparameters$lambda)
  if (isTRUE(x$hyperparameters$enet)) {
    pen_str <- paste0(pen_str, ", lambda1 = ", x$hyperparameters$lambda1)
  }
  cat("Penalty:", pen_str, "\n")

  cat("Constraints: disjoint =", x$hyperparameters$disjoint, "\n")

  invisible(x)
}


#' Summarise an orthoMTL Object
#'
#' Displays a detailed summary of a fitted \code{orthoMTL} model,
#' including hyperparameters, coefficient matrix statistics, and
#' top features by mean absolute coefficient.
#'
#' @param object A fitted model object of class \code{"orthoMTL"}.
#' @param n_top Number of top features to display. Default: 5.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{object}.
#'
#' @method summary orthoMTL
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 100; p <- 10; n_tasks <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' Y <- X %*% matrix(rnorm(p * n_tasks), p) + matrix(rnorm(n * n_tasks), n) * 0.1
#' fit <- orthoMTL(X, Y, lambda = 1e-3)
#' summary(fit)
summary.orthoMTL <- function(object, n_top = 5, ...) {

  # Print compact header
  print.orthoMTL(object)

  cat("---\n")

  # Hyperparameters
  cat("Hyperparameters:\n")
  hp <- object$hyperparameters
  for (nm in names(hp)) {
    cat("  ", nm, "=", hp[[nm]], "\n")
  }

  cat("---\n")

  # Coefficient matrix summary
  B <- object$B
  n_zero <- sum(B == 0)
  n_total <- length(B)
  sparsity <- 100 * n_zero / n_total

  cat("Coefficients:", object$n_features, "x", object$n_tasks, "matrix\n")
  cat("  Sparsity:", format(sparsity, digits = 1, nsmall = 1), "%\n")
  cat("  Range: [", format(min(B), digits = 4), ",",
      format(max(B), digits = 4), "]\n")

  cat("---\n")

  # Top features by mean absolute coefficient
  if (!is.null(object$feature_names)) {
    mean_abs <- sort(apply(abs(B), 1, mean), decreasing = TRUE)
    n_show <- min(n_top, length(mean_abs))
    cat("Top", n_show, "features (mean |coefficient| across tasks):\n")
    top <- mean_abs[seq_len(n_show)]
    for (i in seq_along(top)) {
      cat("  ", names(top)[i], ":", format(top[i], digits = 4), "\n")
    }
  }

  # Survival task thresholds
  if (isTRUE(object$hyperparameters$survival) && !is.null(object$task_names)) {
    cat("---\n")
    cat("Task thresholds:", paste(object$task_names, collapse = ", "), "\n")
  }

  invisible(object)
}
