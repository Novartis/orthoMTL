#' Print Cross-Validation Results
#'
#' Displays a summary of cross-validation results from
#' \code{\link{cv_orthoMTL}}.
#'
#' @param x An object of class \code{"cv_orthoMTL"}.
#' @param n_top Number of top configurations to display. Default: 5.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print cv_orthoMTL
#' @export
#'
#' @examples
#' # See ?cv_orthoMTL for a full example
print.cv_orthoMTL <- function(x, n_top = 5, ...) {

  cat("Cross-validation results for orthoMTL\n")
  cat("Configurations tested:", x$n_configs, "\n")
  cat("Folds:", x$n_folds, "\n")

  cat("---\n")

  cat("Best configuration:\n")
  cat("  lambda   =", x$best$lambda, "\n")
  cat("  alpha    =", x$best$alpha, "\n")
  cat("  stepsize =", x$best$stepsize, "\n")
  cat("  diag_val =", x$best$diag_val, "\n")
  cat("  CV C-index =", format(x$best$cv_score, digits = 4), "\n")

  cat("---\n")

  n_show <- min(n_top, nrow(x$results))
  cat("Top", n_show, "configurations:\n")
  print(x$results[seq_len(n_show), ], row.names = FALSE)

  invisible(x)
}
