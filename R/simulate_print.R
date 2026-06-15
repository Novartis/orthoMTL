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

  mode <- if (!is.null(x$mode)) x$mode else "survival"
  cat("Simulated multi-task ", mode, " data\n", sep = "")
  cat("Patients:", x$n, "| Features:", x$p + 1,
      "(incl. treatment)\n")
  cat("  Continuous:", x$n_continuous,
      "| Binary:", x$p - x$n_continuous,
      "| Treatment: 1\n")
  cat("  Signal features:", x$n_signals,
      "| Null features:", x$p - x$n_signals, "\n")
  cat("Tasks:", paste(x$thresholds, collapse = ", "), "\n")

  if (mode == "survival") {
    event_rate <- round(100 * mean(x$Event), 1)
    cat("Event rate:", event_rate, "%\n")
    cat("Median survival time:",
        format(median(x$SurvTime), digits = 3), "\n")
  } else if (mode == "classification") {
    pos_rate <- round(100 * mean(x$Y > 0), 1)
    cat("Positive-label rate:", pos_rate, "%\n")
  } else {
    cat("Response range: [",
        format(min(x$Y), digits = 3), ", ",
        format(max(x$Y), digits = 3), "]\n", sep = "")
  }

  cat("---\n")
  cat("Signal features and effect types:\n")
  for (feat in x$ground_truth$signal_features) {
    etype <- x$ground_truth$effect_types[feat]
    cat("  ", feat, ":", etype, "\n")
  }

  cat("---\n")
  effect_scale <- if (mode == "survival") "log-hazard" else "coefficient"
  cat("Treatment effect:", x$ground_truth$coefficients["treatment", 1],
      "(", effect_scale, ", constant across tasks)\n")

  invisible(x)
}
