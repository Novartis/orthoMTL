#' Print Bootstrap Inference Results
#'
#' Displays a summary of bootstrap inference results from
#' \code{\link{bootstrap_orthoMTL}}.
#'
#' @param x An object of class \code{"bootstrap_orthoMTL"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print bootstrap_orthoMTL
#' @export
#'
#' @examples
#' # See ?bootstrap_orthoMTL for a full example
print.bootstrap_orthoMTL <- function(x, ...) {

  cat("Bootstrap inference for orthoMTL\n")
  cat("Repeats:", x$n_repeats, "(real) +",
      x$n_repeats, "(null permutation)\n")
  cat("Features:", x$n_features, "| Tasks:", x$n_tasks, "\n")

  cat("---\n")

  cat("Real models -- objective:\n")
  cat("  mean =", format(mean(x$obj_real), digits = 4),
      " range = [", format(min(x$obj_real), digits = 4), ",",
      format(max(x$obj_real), digits = 4), "]\n")

  cat("Null models -- objective:\n")
  cat("  mean =", format(mean(x$obj_null), digits = 4),
      " range = [", format(min(x$obj_null), digits = 4), ",",
      format(max(x$obj_null), digits = 4), "]\n")

  invisible(x)
}
