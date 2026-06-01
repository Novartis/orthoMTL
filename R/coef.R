#' Extract Coefficients from an orthoMTL Model
#'
#' Returns the coefficient matrix from a fitted \code{orthoMTL} model.
#'
#' @param object A fitted model object of class \code{"orthoMTL"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric matrix of regression coefficients with dimensions
#'   \code{n_features x n_tasks}. Row names correspond to feature names
#'   and column names to task names, if available.
#'
#' @method coef orthoMTL
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 100; p <- 10; n_tasks <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' Y <- X %*% matrix(rnorm(p * n_tasks), p) + matrix(rnorm(n * n_tasks), n) * 0.1
#' K <- matrix(1, n_tasks, n_tasks); diag(K) <- 0.5
#' fit <- orthoMTL(X, Y, lambda = 1e-3, K = K)
#' coef(fit)
coef.orthoMTL <- function(object, ...) {
  return(object$B)
}
