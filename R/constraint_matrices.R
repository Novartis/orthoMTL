#' Create Diffusion Constraint Matrix
#'
#' Builds a square constraint matrix \code{K} encoding the prior belief
#' that temporally distant tasks should have more orthogonal coefficients.
#' Off-diagonal entries accumulate weight proportional to the distance
#' between task indices, implementing a diffusion-like pattern.
#'
#' This matrix is passed to \code{\link{orthoMTL}} via the \code{K}
#' argument to control the orthogonality penalty between tasks.
#'
#' @param numTasks An integer specifying the number of tasks.
#' @param diag_val A numeric value for the diagonal of \code{K}.
#'   Default: \code{0.5}. This value is typically overridden during
#'   cross-validation (see \code{\link{cv_orthoMTL}}).
#'
#' @return A numeric square matrix of dimensions
#'   \code{numTasks x numTasks}. Off-diagonal entry \code{K[i,j]}
#'   is larger when tasks \code{i} and \code{j} are further apart.
#'   Diagonal entries are set to \code{diag_val}.
#'
#' @details
#' The construction rule: for each distance level \code{h} from
#' \code{1} to \code{numTasks}, add \code{1} to all entries where
#' \code{|row - col| > h}. This produces a matrix where nearby tasks
#' (adjacent thresholds) share more support, while distant tasks
#' are pushed toward orthogonality.
#'
#' For survival analysis with time thresholds, this encodes the
#' assumption that the set of predictive features changes gradually
#' over time rather than abruptly.
#'
#' @seealso \code{\link{cv_orthoMTL}}
#'
#' @export
#'
#' @examples
#' # 5-task constraint matrix
#' K <- create_constraint_matrix(5)
#' K
#'
#' # Override diagonal for a specific penalty balance
#' K <- create_constraint_matrix(7, diag_val = 6)
#' K
create_constraint_matrix <- function(numTasks, diag_val = 0.5) {

  K <- matrix(0, nrow = numTasks, ncol = numTasks)
  delta <- row(K) - col(K)

  # Accumulate off-diagonal weights by distance
  for (high in seq_len(numTasks)) {
    K[abs(delta) > high] <- K[abs(delta) > high] + 1
  }

  # Set diagonal
  diag(K) <- diag_val

  return(K)
}
