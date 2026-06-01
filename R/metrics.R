#' Concordance Index for Multi-Task Predictions
#'
#' Computes a concordance index (C-index) adapted for the multi-task
#' survival framework. Predictions across tasks are aggregated by
#' row-sum to produce a single score per patient, then concordance
#' is evaluated over pairs where at least one member has a fully
#' observed outcome.
#'
#' @param true.label.mat A numeric matrix of true response labels with
#'   dimensions \code{n x numTasks}. May contain \code{NA} for censored
#'   observations (as produced by \code{\link{create_longitudinal_labels}}).
#' @param pred.label.mat A numeric matrix of predicted response values
#'   with dimensions \code{n x numTasks} (as produced by
#'   \code{\link{predict.orthoMTL}}).
#'
#' @return A numeric value between 0 and 1 (higher is better).
#'   A value of 0.5 indicates random concordance.
#'
#' @details
#' The effective survival time for each patient is derived as the
#' highest task index where the true label is positive (i.e., the
#' last threshold at which the patient was progression-free).
#'
#' A patient is considered "uncensored" only if all task labels are
#' non-\code{NA}. Concordant pairs require both a correct ordering
#' of effective survival times and a matching ordering of predicted
#' scores.
#'
#' @section Known limitations:
#' The following limitations are documented and flagged for future
#' investigation:
#' \itemize{
#'   \item No credit for tied predictions or tied survival times
#'   \item Strict censoring: uncensored requires all tasks observed
#'   \item Equal weighting of all tasks in the row-sum aggregation
#' }
#'
#' @seealso \code{\link{predict.orthoMTL}} for generating the
#'   prediction matrix.
#'
#' @export
#'
#' @examples
#' # Simulate a small multi-task prediction scenario
#' set.seed(42)
#' n <- 30; n_tasks <- 4
#' SurvTime <- rexp(n, rate = 0.1)
#' Event <- rbinom(n, 1, 0.7)
#' thresholds <- c(4, 6, 10, 15)
#'
#' Y <- create_longitudinal_labels(SurvTime, Event, thresholds)
#'
#' # Simulate imperfect predictions (add noise to true labels)
#' pred <- Y
#' pred[is.na(pred)] <- 0.5
#' pred <- pred + matrix(rnorm(n * n_tasks, sd = 0.3), n, n_tasks)
#'
#' cindex_mtl(Y, pred)
cindex_mtl <- function(true.label.mat, pred.label.mat) {

  # TODO(v1.1): Row-sum aggregation weights all tasks equally regardless
  #   of interval width (e.g., threshold gap 2->4 months counts the same
  #   as 15->22 months). Investigate whether weighted aggregation by
  #   interval width improves discrimination.
  pred.vec <- apply(pred.label.mat, 1, sum)

  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  sum4 <- 0

  for (i in 1:(nrow(true.label.mat) - 1)) {
    for (j in (i + 1):nrow(true.label.mat)) {

      # Effective survival time: last threshold where label > 0
      if (length(which(true.label.mat[i, ] > 0)) > 0) {
        stime1 <- max(which(true.label.mat[i, ] > 0))
      } else {
        stime1 <- 0
      }

      if (length(which(true.label.mat[j, ] > 0)) > 0) {
        stime2 <- max(which(true.label.mat[j, ] > 0))
      } else {
        stime2 <- 0
      }

      pred1 <- pred.vec[i]
      pred2 <- pred.vec[j]

      # TODO(v1.1): Strict censoring definition — a patient is considered
      #   uncensored only if ALL task labels are non-NA. Patients censored
      #   at later thresholds never contribute as the event member of a
      #   concordant pair. This may bias the metric toward patients with
      #   early events. Investigate per-task censoring alternatives.
      status1 <- as.numeric(all(!is.na(true.label.mat[i, ])))
      status2 <- as.numeric(all(!is.na(true.label.mat[j, ])))

      # TODO(v1.1): No tie handling — concordance only counts strict
      #   inequalities. Standard implementations (e.g., survival::concordance)
      #   give 0.5 credit for ties. With binary features producing many
      #   tied predictions, this may systematically undercount concordance.
      #   Investigate impact on reported C-index values.
      if (stime1 < stime2 & pred1 < pred2 & status1 == 1) sum1 <- sum1 + 1
      if (stime2 < stime1 & pred2 < pred1 & status2 == 1) sum2 <- sum2 + 1
      if (stime1 < stime2 & status1 == 1) sum3 <- sum3 + 1
      if (stime2 < stime1 & status2 == 1) sum4 <- sum4 + 1
    }
  }

  return((sum1 + sum2) / (sum3 + sum4))
}
