#' Print Simulated Multi-Task Data
#'
#' Displays a summary of a simulated dataset from
#' \code{\link{simulate_mtl}}.
#'
#' @param x An object of class \code{"simulated_mtl"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print simulated_mtl
#' @export
#'
#' @examples
#' sim <- simulate_mtl(n = 100, p = 20, n_signals = 4)
#' print(sim)
print.simulated_mtl <- function(x, ...) {

  cat("Simulated multi-task survival data\n")
  cat("Patients:", x$n, "| Features:", x$p + 1,
      "(incl. treatment)\n")
  cat("  Continuous:", x$n_continuous,
      "| Binary:", x$p - x$n_continuous,
      "| Treatment: 1\n")
  cat("  Signal features:", x$n_signals,
      "| Null features:", x$p - x$n_signals, "\n")
  cat("Thresholds:", paste(x$thresholds, collapse = ", "), "\n")

  event_rate <- round(100 * mean(x$Event), 1)
  cat("Event rate:", event_rate, "%\n")

  cat("Median survival time:",
      format(median(x$SurvTime), digits = 3), "\n")

  cat("---\n")
  cat("Signal features and effect types:\n")
  for (feat in x$ground_truth$signal_features) {
    etype <- x$ground_truth$effect_types[feat]
    cat("  ", feat, ":", etype, "\n")
  }

  cat("---\n")
  cat("Treatment effect:", x$ground_truth$coefficients["treatment", 1],
      "(log-hazard, constant across thresholds)\n")

  invisible(x)
}
