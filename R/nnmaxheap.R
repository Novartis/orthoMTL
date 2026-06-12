#' Project a Vector onto Non-Negative Non-Increasing Space
#'
#' Applies an isotonic regression-style projection to enforce that the
#' output vector is non-negative and non-increasing. Used internally
#' by \code{\link{predict.orthoMTL}} to ensure survival predictions
#' are monotonically decreasing across time thresholds.
#'
#' @param m A numeric vector to project.
#'
#' @return A numeric vector of the same length as \code{m}, projected
#'   onto the non-negative non-increasing constraint space.
#'
#' @details
#' This function implements a pool-adjacent-violators style algorithm.
#' It replaces the external \code{Iso} package dependency used in the
#' predecessor \code{orthopen} package.
#'
#' @export
#'
#' @examples
#' # Project a vector onto the non-negative non-increasing space
#' nnmaxheap_C(c(3, 1, 2, -1))
#'
#' # Already valid input: returned unchanged
#' nnmaxheap_C(c(5, 3, 3, 1))
nnmaxheap_C <- function(m) {
  n <- length(m)

  if (n < 1) {
    stop("n=", n, " should be an integer over 1!")
  }

  # handle length-1 case
  if (n == 1) {
    return(max(m, 0))   # non-negative projection of single value
  }

  x <- rep(0, n)
  location <- rep(0, n)

  i <- n
  x[i] <- m[i]
  location[i] <- i

  # Process remaining elements bottom-up
  for (i in (n - 1):1) {
    if (m[i] > x[i + 1]) {
      x[i] <- m[i]
      location[i] <- i
    } else {
      # Merge with the first group
      num <- location[i + 1] - i
      x[i] <- (m[i] + x[i + 1] * num) / (num + 1)
      location[i] <- location[i + 1]
      j <- location[i + 1] + 1
      while (TRUE) {
        if (j > n) break
        if (x[i] <= x[j]) {
          num <- location[j] - j + 1
          x[i] <- (x[i] * (j - i) + x[j] * num) / (location[j] - i + 1)
          location[i] <- location[j]
          j <- location[j] + 1
        } else {
          break
        }
      }
    }
  }

  # Compute the solution using mean and location
  i <- 1
  while (i < n) {
    if (x[i] > 0) {
      for (j in (i + 1):location[i]) {
        x[j] <- x[i]
      }
      i <- location[i] + 1
    } else {
      for (j in i:n) {
        x[j] <- 0
      }
      break
    }
  }

  return(x)
}
