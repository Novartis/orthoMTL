#' @title Projection onto Disjoint Support Feasible Set
#' @description Projects weight and constraint matrices ...
#' @param w Numeric matrix of weights used in the loss function.
#' @param v Numeric matrix of constraints used for projection.
#' @return A list with projected matrices \code{w} and \code{v}.
#' @keywords internal
proj_disjoint <- function(w, v) {
  w_tmp = w
  v_tmp = v
  idx1 = v <= 0
  v_tmp[idx1] = 0
  w_tmp[idx1] = 0
  abs_w = abs(w)
  idx2 = abs_w < v
  idx3 = !(idx1 | idx2)
  tmp = 0.5 * (v[idx3] + abs_w[idx3])
  v_tmp[idx3] = tmp
  w_tmp[idx3] = sign(w[idx3]) * tmp
  return(list(w = w_tmp, v = v_tmp))
}
