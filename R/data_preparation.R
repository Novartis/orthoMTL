#' Conversion of Survival data into Longitudinal data
#'
#' @param SurvTime a vector of n observed times
#' @param Event a binary vector encoding censoring of the data for n patients.
#'   0 = censored, 1 = event was observed.
#' @param thresholds a numeric vector of length numTasks representing the
#'   considered cut-off values
#'
#' @return Y a matrix of size n x numTasks representing the longitudinal
#'   censored data
#' @export
create_longitudinal_labels <- function(SurvTime, Event, thresholds = c(4, 6)) {

  Y <- matrix(1, nrow = length(SurvTime), ncol = length(thresholds))

  for (i in seq_along(thresholds)) {
    Y[which(SurvTime < thresholds[i] & Event == 1), i] <- 0
  }

  for (i in seq_along(thresholds)) {
    Y[which(SurvTime < thresholds[i] & Event == 0), i] <- NA
  }

  return(Y)
}


#' Indicator matrix for censored data
#'
#' @param Y a matrix of size n x numTasks representing the longitudinal
#'   censored data
#'
#' @return W a matrix of size n x numTasks where 1 = observed, 0 = censored
#' @export
create_indicator_matrix <- function(Y) {
  W <- matrix(1, nrow = nrow(Y), ncol = ncol(Y))
  W[which(is.na(Y), arr.ind = TRUE)] <- 0
  return(W)
}
