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


# ---------------------------------------------------------------------
# General-purpose metrics for non-survival modes (GEN-01)
#
# These complement the survival-specific cindex_mtl for the regression
# and classification use cases of orthoMTL. All operate cell-wise on the
# n x numTasks true/predicted matrices, pooling over every non-NA cell
# (so a single NA-masked or fully observed matrix both work). They never
# touch cindex_mtl, which stays faithful to the MTLSA reference.
# ---------------------------------------------------------------------

# Internal: flatten true/pred matrices to aligned non-NA vectors.
.mtl_pool <- function(true.mat, pred.mat) {
  true.mat <- as.matrix(true.mat)
  pred.mat <- as.matrix(pred.mat)
  if (!all(dim(true.mat) == dim(pred.mat))) {
    stop("true and predicted matrices must have identical dimensions.",
         call. = FALSE)
  }
  t_vec <- as.numeric(true.mat)
  p_vec <- as.numeric(pred.mat)
  keep <- !is.na(t_vec) & !is.na(p_vec)
  if (!any(keep)) {
    stop("No non-NA cells to evaluate.", call. = FALSE)
  }
  list(true = t_vec[keep], pred = p_vec[keep])
}

#' Root Mean Squared Error for Multi-Task Regression
#'
#' Pooled RMSE over every non-\code{NA} cell of the true/predicted
#' matrices, for \code{orthoMTL} fits in regression mode
#' (\code{logistic = FALSE}, \code{survival = FALSE}).
#'
#' @param true.label.mat A numeric matrix of true responses
#'   (\code{n x numTasks}). \code{NA} cells are ignored.
#' @param pred.label.mat A numeric matrix of predictions of the same
#'   dimensions (as produced by \code{\link{predict.orthoMTL}}).
#'
#' @return A single non-negative numeric value (lower is better).
#'
#' @seealso \code{\link{r2_mtl}}, \code{\link{cindex_mtl}}
#' @export
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(30), 10, 3)
#' P <- Y + matrix(rnorm(30, sd = 0.1), 10, 3)
#' rmse_mtl(Y, P)
rmse_mtl <- function(true.label.mat, pred.label.mat) {
  v <- .mtl_pool(true.label.mat, pred.label.mat)
  sqrt(mean((v$true - v$pred)^2))
}

#' Coefficient of Determination (R-squared) for Multi-Task Regression
#'
#' Pooled \eqn{R^2 = 1 - SS_{res} / SS_{tot}} over every non-\code{NA}
#' cell, for \code{orthoMTL} fits in regression mode. \eqn{SS_{tot}} uses
#' the global mean of the pooled true values.
#'
#' @inheritParams rmse_mtl
#'
#' @return A single numeric value (higher is better; 1 = perfect).
#'   Can be negative when predictions are worse than the mean.
#'
#' @seealso \code{\link{rmse_mtl}}, \code{\link{cindex_mtl}}
#' @export
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(30), 10, 3)
#' P <- Y + matrix(rnorm(30, sd = 0.1), 10, 3)
#' r2_mtl(Y, P)
r2_mtl <- function(true.label.mat, pred.label.mat) {
  v <- .mtl_pool(true.label.mat, pred.label.mat)
  ss_res <- sum((v$true - v$pred)^2)
  ss_tot <- sum((v$true - mean(v$true))^2)
  if (ss_tot == 0) {
    stop("Cannot compute R-squared: true values have zero variance.",
         call. = FALSE)
  }
  1 - ss_res / ss_tot
}

#' Classification Accuracy for Multi-Task Predictions
#'
#' Pooled fraction of correctly classified non-\code{NA} cells, for
#' \code{orthoMTL} fits in classification mode (\code{logistic = TRUE}).
#' A cell is positive when its true label is \code{> 0} and predicted
#' positive when its score is \code{> 0}. This works for both the
#' \eqn{\{-1, +1\}} encoding the solver optimises and a \eqn{\{0, 1\}}
#' encoding, and for raw \code{"link"} scores or \code{"response"}
#' probabilities (threshold 0.5 corresponds to a link of 0... see note).
#'
#' @inheritParams rmse_mtl
#' @param threshold Decision threshold applied to \code{pred.label.mat}.
#'   Default \code{0} matches raw link scores; pass \code{0.5} for
#'   probabilities from \code{predict(..., type = "response")}.
#'
#' @return A single numeric value in \code{[0, 1]} (higher is better).
#'
#' @seealso \code{\link{auc_mtl}}, \code{\link{cindex_mtl}}
#' @export
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(sample(c(-1, 1), 30, replace = TRUE), 10, 3)
#' P <- Y * abs(matrix(rnorm(30), 10, 3))   # mostly correct signs
#' accuracy_mtl(Y, P)
accuracy_mtl <- function(true.label.mat, pred.label.mat, threshold = 0) {
  v <- .mtl_pool(true.label.mat, pred.label.mat)
  true_pos <- v$true > 0
  pred_pos <- v$pred > threshold
  mean(true_pos == pred_pos)
}

#' Area Under the ROC Curve for Multi-Task Predictions
#'
#' Pooled AUC over every non-\code{NA} cell, computed via the
#' Mann-Whitney U statistic (with 0.5 credit for ties), for
#' \code{orthoMTL} fits in classification mode. The positive class is
#' defined by true label \code{> 0}. AUC is invariant to monotone
#' transforms, so raw \code{"link"} scores and \code{"response"}
#' probabilities give identical results.
#'
#' @inheritParams rmse_mtl
#'
#' @return A single numeric value in \code{[0, 1]} (higher is better;
#'   0.5 = random).
#'
#' @seealso \code{\link{accuracy_mtl}}, \code{\link{cindex_mtl}}
#' @export
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(sample(c(-1, 1), 30, replace = TRUE), 10, 3)
#' P <- Y + matrix(rnorm(30), 10, 3)
#' auc_mtl(Y, P)
auc_mtl <- function(true.label.mat, pred.label.mat) {
  v <- .mtl_pool(true.label.mat, pred.label.mat)
  pos <- v$pred[v$true > 0]
  neg <- v$pred[v$true <= 0]
  n_pos <- length(pos)
  n_neg <- length(neg)
  if (n_pos == 0 || n_neg == 0) {
    stop("AUC requires both positive (label > 0) and negative cells.",
         call. = FALSE)
  }
  # Mann-Whitney U via rank-sum, with 0.5 credit for ties.
  r <- rank(c(pos, neg))
  auc <- (sum(r[seq_len(n_pos)]) - n_pos * (n_pos + 1) / 2) /
    (n_pos * n_neg)
  auc
}
