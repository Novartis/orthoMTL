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


#' Check Whether a Constraint Matrix Gives a Convex Penalty
#'
#' Reports whether the orthogonality penalty
#' \eqn{\Omega_K(W)^2 = \sum_{s,t} K_{st} |W_s^\top W_t|} defined by a
#' constraint matrix \code{K} is convex, and how far its diagonal is from
#' the convexity threshold.
#'
#' @param K A square, symmetric, non-negative numeric matrix (one row and
#'   column per task), as passed to \code{\link{orthoMTL}}.
#' @param tol Numerical tolerance used for the symmetry check and on the
#'   smallest eigenvalue, relative to \code{max(1, max(abs(K)))}.
#'   Default: \code{1e-10}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{convex}}{Logical: \code{TRUE} if the penalty is convex.}
#'     \item{\code{min_eigenvalue}}{Smallest eigenvalue of the companion
#'       matrix (diagonal of \code{K} kept, off-diagonal entries negated).
#'       The penalty is convex if and only if this is \eqn{\geq 0}.}
#'     \item{\code{lambda_max_offdiag}}{Largest eigenvalue of the
#'       off-diagonal part of \code{K} (diagonal set to 0). This is the
#'       smallest constant diagonal that makes \code{K} convex.}
#'     \item{\code{rho}}{\code{min(diag(K)) / lambda_max_offdiag}. Values
#'       \code{< 1} are non-convex, \code{1} is the atomic-norm limit and
#'       \code{> 1} is strictly convex. \code{Inf} when \code{K} is
#'       diagonal.}
#'   }
#'
#' @details
#' Write \eqn{K = O + D}, with \eqn{O} the off-diagonal part and \eqn{D} the
#' diagonal. For non-negative \code{K}, the penalty is convex if and only
#' if the companion matrix \eqn{D - O} is positive semi-definite (Vervier et
#' al., 2014). With a constant diagonal \eqn{D = dI}, this is
#' \eqn{d \geq \lambda_{\max}(O)}, i.e. \code{rho >= 1}.
#'
#' A non-convex \code{K} (\code{rho < 1}) applies a stronger orthogonality
#' pressure. That can help when the orthogonal task pairs encoded in
#' \code{K} are the true ones, but it is costly when they are not. For
#' example, \code{create_constraint_matrix()} with its default diagonal is
#' strongly non-convex.
#'
#' @references Vervier, K., Mahé, P., d'Aspremont, A., Veyrieras, J.-B., &
#'   Vert, J.-P. (2014). On Learning Matrices with Orthogonal Columns or
#'   Disjoint Supports. \emph{ECML-PKDD 2014}.
#'
#' @seealso \code{\link{create_constraint_matrix}}, \code{\link{orthoMTL}}
#'
#' @export
#'
#' @examples
#' # Identity (per-task ridge): convex, no coupling between tasks
#' validate_K(diag(4))$convex
#'
#' # The diffusion matrix with its default diagonal is non-convex
#' K <- create_constraint_matrix(5)
#' chk <- validate_K(K)
#' chk$convex
#' chk$rho
#'
#' # Smallest constant diagonal that makes it convex
#' K_convex <- create_constraint_matrix(5, diag_val = chk$lambda_max_offdiag)
#' validate_K(K_convex)$convex
validate_K <- function(K, tol = 1e-10) {

  if (!is.matrix(K) || !is.numeric(K) || nrow(K) != ncol(K)) {
    stop("'K' must be a square numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(K))) {
    stop("'K' must not contain NA, NaN or infinite values.", call. = FALSE)
  }
  scale <- max(1, max(abs(K)))
  if (max(abs(K - t(K))) > tol * scale) {
    stop("'K' must be symmetric.", call. = FALSE)
  }
  if (any(K < 0)) {
    stop("'K' must be non-negative.", call. = FALSE)
  }

  d <- diag(K)
  O <- K
  diag(O) <- 0

  lambda_max_off <- max(eigen(O, symmetric = TRUE, only.values = TRUE)$values)
  companion <- diag(d, nrow = length(d)) - O
  min_eig <- min(eigen(companion, symmetric = TRUE, only.values = TRUE)$values)

  rho <- if (lambda_max_off > tol * scale) min(d) / lambda_max_off else Inf

  list(
    convex             = min_eig >= -tol * scale,
    min_eigenvalue     = min_eig,
    lambda_max_offdiag = lambda_max_off,
    rho                = rho
  )
}
