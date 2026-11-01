#' Max-cut semidefinite program for the Omega_O penalty (internal)
#'
#' Solves \eqn{\max \{ \langle G, Z \rangle : Z \succeq 0,\ \mathrm{diag}(Z) = 1 \}}
#' for a small positive semidefinite \eqn{T \times T} matrix \code{G} by the
#' Burer-Monteiro "mixing method" (block-coordinate ascent on unit vectors
#' \eqn{Z = V^\top V}, rank \eqn{T}), and certifies the result with a dual
#' feasible point.
#'
#' For a coefficient matrix \eqn{W} with \eqn{p \ge T} rows, the optimal value
#' with \eqn{G = W^\top W} is \eqn{\Omega_O(W)^2}, the squared atomic norm of
#' matrices with orthogonal columns (equivalently the squared 2-summing norm of
#' \eqn{W : \ell_\infty^T \to \ell_2^p}; Tropp 2009, Sec. 3.3), and
#' \eqn{2 W Z^\star} is a subgradient of \eqn{\Omega_O(W)^2} at \eqn{W}.
#' For \eqn{p < T} the value is a convex lower bound of \eqn{\Omega_O(W)^2}.
#'
#' @param G A symmetric positive semidefinite matrix.
#' @param V0 Optional warm start, a \eqn{T \times T} matrix with unit columns;
#'   it is re-inflated by \code{inflate} \eqn{\times I} before use so that a
#'   collapsed (low-rank) warm start cannot trap the iteration below the SDP
#'   optimum.
#' @param inflate Size of that re-inflation. Inside the solver loop a tiny
#'   value suffices: the escape from a collapsed start then continues across
#'   successive warm-started calls.
#'   Default: the identity (deterministic; no random numbers are drawn).
#' @param tol Relative tolerance on the objective between two sweeps.
#' @param max_sweeps Maximum number of coordinate sweeps.
#' @param certify Compute the dual certificate (one eigendecomposition)?
#'   \code{FALSE} inside the solver loop, where only the value and \code{Z}
#'   are needed; \code{gap} is then \code{NA}.
#' @return A list with \code{value} (the primal value
#'   \eqn{\langle G, Z \rangle}), \code{Z}, \code{V} (for warm starts) and
#'   \code{gap} (an upper bound on the optimality gap, from the dual point
#'   \eqn{\Lambda = \mathrm{diag}(GZ) + \mu I} with \eqn{\mu} the smallest
#'   shift making \eqn{\Lambda - G \succeq 0}).
#' @references Joel A. Tropp (2009). Column subset selection, matrix
#'   factorization, and eigenvalue optimization. SODA. arXiv:0806.4404.
#' @keywords internal
#' @noRd
.maxcut_sdp <- function(G, V0 = NULL, tol = 1e-10, max_sweeps = 500,
                        certify = TRUE, inflate = 1e-3) {
  n <- nrow(G)
  if (is.null(V0)) {
    V <- diag(1, n)
  } else {
    # The mixing update keeps V inside its column span, so a warm start that
    # has collapsed to low rank (rank 1 is optimal for a balanced G) could never
    # reach a higher-rank optimum for a later, frustrated G. Re-inflate it
    # deterministically (no random numbers) and renormalise the columns.
    V <- V0 + inflate * diag(1, n)
    V <- sweep(V, 2, sqrt(colSums(V^2)), "/")
  }
  val_old <- sum(G * crossprod(V))
  for (sweep in seq_len(max_sweeps)) {
    for (i in seq_len(n)) {
      # best unit vector for column i given the others: maximise <v_i, sum_j G_ij v_j>
      g <- V[, -i, drop = FALSE] %*% G[-i, i]
      ng <- sqrt(sum(g^2))
      if (ng > 0) V[, i] <- g / ng
    }
    val <- sum(G * crossprod(V))
    if (abs(val - val_old) <= tol * max(1, abs(val))) break
    val_old <- val
  }
  Z <- crossprod(V)
  if (!certify) return(list(value = val, Z = Z, V = V, gap = NA_real_))
  # dual certificate: Lambda = diag(G Z) + mu I with Lambda - G psd
  lam <- rowSums(G * Z)
  mu <- max(0, -min(eigen(diag(lam, n) - G, symmetric = TRUE,
                          only.values = TRUE)$values))
  gap <- sum(lam) + n * mu - val
  # second dual point, c * diag(lam) with c = lambda_max(D^-1/2 G D^-1/2):
  # often tighter than the additive shift; keep the better bound
  if (all(lam > 0)) {
    d <- 1 / sqrt(lam)
    cmul <- max(eigen(outer(d, d) * G, symmetric = TRUE, only.values = TRUE)$values)
    gap <- min(gap, cmul * sum(lam) - val)
  }
  list(value = val, Z = Z, V = V, gap = gap)
}
