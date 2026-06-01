#' Convert Survival Data to Longitudinal Binary Labels
#'
#' Transforms time-to-event survival data into a matrix of binary labels
#' at user-defined time thresholds. This is the core data transformation
#' that enables survival analysis within the multi-task learning framework.
#'
#' The encoding logic for each patient at each threshold is:
#' \itemize{
#'   \item \code{1} — patient is progression-free at this threshold
#'     (survival time exceeds threshold, regardless of event status)
#'   \item \code{0} — patient experienced an event before this threshold
#'     (survival time < threshold AND event observed)
#'   \item \code{NA} — patient was censored before this threshold
#'     (survival time < threshold AND no event observed); label is unknown
#' }
#'
#' @param SurvTime A numeric vector of length \code{n} containing
#'   observed survival times. Must be non-negative.
#' @param Event A numeric vector of length \code{n} encoding censoring
#'   status. \code{1} = event observed, \code{0} = censored.
#' @param thresholds A numeric vector of length \code{numTasks}
#'   representing the time thresholds at which to evaluate survival
#'   status. Default: \code{c(4, 6)}.
#'
#' @return A numeric matrix of dimensions \code{n x numTasks}. Column
#'   names are set to the threshold values. Contains \code{1}, \code{0},
#'   and \code{NA} values as described above.
#'
#' @seealso \code{\link{create_indicator_matrix}} to convert \code{NA}
#'   values into a binary censoring indicator matrix.
#'
#' @export
#'
#' @examples
#' # Simulate 10 patients
#' set.seed(42)
#' SurvTime <- rexp(10, rate = 0.1)
#' Event <- rbinom(10, 1, 0.7)
#' thresholds <- c(4, 6, 10, 15)
#'
#' Y <- create_longitudinal_labels(SurvTime, Event, thresholds)
#' Y   # 1 = progression-free, 0 = event, NA = censored
create_longitudinal_labels <- function(SurvTime, Event, thresholds = c(4, 6)) {

  if (any(SurvTime < 0)) {
    stop("SurvTime vector should be positive values only", call. = FALSE)
  }
  if (any(Event > 1)) {
    stop("Event vector is encoded as binary values", call. = FALSE)
  }

  Y <- matrix(1, nrow = length(SurvTime), ncol = length(thresholds))

  # Observed events before each threshold
  for (i in seq_along(thresholds)) {
    Y[which(SurvTime < thresholds[i] & Event == 1), i] <- 0
  }

  # Censored observations before each threshold
  for (i in seq_along(thresholds)) {
    Y[which(SurvTime < thresholds[i] & Event == 0), i] <- NA
  }

  colnames(Y) <- thresholds
  return(Y)
}


#' Create Censoring Indicator Matrix
#'
#' Converts \code{NA} entries in a longitudinal label matrix (as produced
#' by \code{\link{create_longitudinal_labels}}) into a binary indicator
#' matrix. This indicator is used by \code{\link{orthoMTL}} to mask
#' censored observations in the loss computation.
#'
#' @param Y A numeric matrix of dimensions \code{n x numTasks}, typically
#'   produced by \code{\link{create_longitudinal_labels}}. May contain
#'   \code{NA} values for censored observations.
#'
#' @return A numeric matrix of the same dimensions as \code{Y}.
#'   \code{1} = label is observed (known), \code{0} = label is censored
#'   (unknown).
#'
#' @seealso \code{\link{create_longitudinal_labels}}
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' SurvTime <- rexp(10, rate = 0.1)
#' Event <- rbinom(10, 1, 0.7)
#' Y <- create_longitudinal_labels(SurvTime, Event, c(4, 6, 10))
#' W <- create_indicator_matrix(Y)
#' W   # 1 = observed, 0 = censored
create_indicator_matrix <- function(Y) {
  W <- matrix(1, nrow = nrow(Y), ncol = ncol(Y))
  W[which(is.na(Y), arr.ind = TRUE)] <- 0
  return(W)
}
