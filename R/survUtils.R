#' C-index for non-Cox model: compare predicted values for non-censored pairs
#'
#' @param true.label.mat a matrix of true response variables with dimensions
#'   n x numTasks. NAs can be used for censored data.
#' @param pred.label.mat a matrix of predicted response variables with the
#'   same dimensions as true.label.mat.
#'
#' @return c-index between 0 and 1 (higher is better)
#' @export
cindex_mtl <- function(true.label.mat, pred.label.mat) {

  pred.vec <- apply(pred.label.mat, 1, sum)
  sum1 <- sum2 <- sum3 <- sum4 <- 0

  for (i in 1:(nrow(true.label.mat) - 1)) {
    for (j in (i + 1):nrow(true.label.mat)) {
      stime1 <- if (any(true.label.mat[i, ] > 0, na.rm = TRUE)) {
        max(which(true.label.mat[i, ] > 0))
      } else 0

      stime2 <- if (any(true.label.mat[j, ] > 0, na.rm = TRUE)) {
        max(which(true.label.mat[j, ] > 0))
      } else 0

      pred1 <- pred.vec[i]
      pred2 <- pred.vec[j]
      status1 <- as.numeric(all(!is.na(true.label.mat[i, ])))
      status2 <- as.numeric(all(!is.na(true.label.mat[j, ])))

      if (stime1 < stime2 & pred1 < pred2 & status1 == 1) sum1 <- sum1 + 1
      if (stime2 < stime1 & pred2 < pred1 & status2 == 1) sum2 <- sum2 + 1
      if (stime1 < stime2 & status1 == 1) sum3 <- sum3 + 1
      if (stime2 < stime1 & status2 == 1) sum4 <- sum4 + 1
    }
  }
  return((sum1 + sum2) / (sum3 + sum4))
}
