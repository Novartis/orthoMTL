#' Predict from an orthoMTL Model
#'
#' Generate predictions from a fitted \code{orthoMTL} model for new
#' observations. In survival mode, predictions are projected onto the
#' non-negative non-increasing space to ensure monotonicity across tasks.
#'
#' @param object A fitted model object of class \code{"orthoMTL"}.
#' @param newdata A numeric matrix of new observations with dimensions
#'   \code{n_new x p}. Column names are used for feature alignment if
#'   available in the fitted object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric matrix of predictions with dimensions
#'   \code{n_new x n_tasks}. Column names correspond to task names
#'   if available.
#'
#' @details
#' If the fitted object contains \code{feature_names} (i.e., the training
#' matrix \code{X} had column names), \code{newdata} columns are aligned
#' to match. Missing features cause an error. Extra features trigger a
#' warning and are dropped.
#'
#' In survival mode (\code{object$hyperparameters$survival == TRUE}),
#' each row of the raw prediction matrix is projected via
#' \code{nnmaxheap_C()} to enforce non-negative, non-increasing values
#' across tasks (time thresholds).
#'
#' @references
#' Vervier, K., Mahe, P., d'Aspremont, A., Veyrieras, J.-B., and
#' Vert, J.-P. (2014). On Learning Matrices with Orthogonal Columns
#' or Disjoint Supports. \emph{ECML-PKDD 2014}.
#' \url{https://hal.science/hal-00985654}
#'
#' @method predict orthoMTL
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 100; p <- 10; n_tasks <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' W_true <- qr.Q(qr(matrix(rnorm(p * n_tasks), p, n_tasks)))
#' Y <- X %*% W_true + matrix(rnorm(n * n_tasks), n) * 0.1
#' K <- matrix(1, n_tasks, n_tasks); diag(K) <- 0.5
#' fit <- orthoMTL(X, Y, lambda = 1e-3, K = K, disjoint = FALSE)
#'
#' X_new <- matrix(rnorm(20 * p), 20, p)
#' colnames(X_new) <- paste0("V", seq_len(p))
#' preds <- predict(fit, newdata = X_new)
#' dim(preds)
predict.orthoMTL <- function(object, newdata, ...) {

  # --- Column alignment ---
  if (!is.null(object$feature_names)) {

    if (is.null(colnames(newdata))) {
      stop("'newdata' must have column names when the fitted model ",
           "contains feature_names.", call. = FALSE)
    }

    expected <- object$feature_names
    provided <- colnames(newdata)

    # Error on missing features
    missing_feats <- setdiff(expected, provided)
    if (length(missing_feats) > 0) {
      stop("The following features required by the model are missing ",
           "from 'newdata': ",
           paste(missing_feats, collapse = ", "), call. = FALSE)
    }

    # Warn on extra features
    extra_feats <- setdiff(provided, expected)
    if (length(extra_feats) > 0) {
      warning("The following features in 'newdata' are not used by the ",
              "model and will be ignored: ",
              paste(extra_feats, collapse = ", "), call. = FALSE)
    }

    # Align and subset
    newdata <- newdata[, expected, drop = FALSE]

  } else {
    # No feature names — fall back to dimension check
    if (ncol(newdata) != object$n_features) {
      stop("'newdata' has ", ncol(newdata), " columns but the model ",
           "expects ", object$n_features, ".", call. = FALSE)
    }
  }

  # --- Prediction ---
  B <- object$B
  Mt <- newdata %*% B

  # --- Survival monotonicity projection ---
  # TODO(v1.1): Monotonicity projection is applied here in predict() but
  #   not during training loss computation. This asymmetry means the model
  #   is trained on raw scores but evaluated on projected scores.
  #   Investigate impact on SOLAR-1 results and simulated vignette.
  if (isTRUE(object$hyperparameters$survival)) {
    Mt <- t(apply(Mt, 1, nnmaxheap_C))
  }

  # --- Label output ---
  if (!is.null(object$task_names)) {
    colnames(Mt) <- object$task_names
  }
  if (!is.null(rownames(newdata))) {
    rownames(Mt) <- rownames(newdata)
  }

  return(Mt)
}
