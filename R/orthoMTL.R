#' Multi-task learning with orthogonal constraints
#'
#' This function solves a multi-task problem where relationships between tasks can be complex
#'
#' @param X a matrix of predictor variables with dimensions n x p
#' @param Y a matrix of response variables with dimensions n x numTasks, where numTasks is the number of response variables. NAs can be used for censored data.
#' @param lambda the regularization parameter for the OrthoPen penalty, default is 1
#' @param step_size the step size for updating the regression coefficients in gradient descent, default is 0.1
#' @param verbose the level of verbosity, default is 0 (no messages)
#' @param stop_no_improve the number of iterations without improvement in the objective function to trigger convergence, default is 100
#' @param max_iter the maximum number of iterations, default is 1e+06
#' @param W_0 a matrix of initial values for the regression coefficients, default is NULL, as not provided and will be randomly attributed
#' @param seed a random seed for reproducibility, default is 42
#' @param K a constraint matrix of weights to adjust the OrthoPen penalty, default is an identity matrix with dimensions numTasks x numTasks
#' @param disjoint a logical value indicating whether the response variables should have disjoint supports, default is FALSE
#' @param logistic a logical value indicating whether logistic regression should be used instead of linear regression, default is FALSE
#' @param alpha the elastic-net mixing parameter in \eqn{[0, 1]}, default is 0.
#'   The penalty is \eqn{\lambda [(1 - \alpha)/2\, \Omega_K(W)^2 + \alpha \|W\|_1]};
#'   \code{alpha = 0} gives the pure orthogonality penalty and \code{alpha = 1}
#'   gives a pure Lasso.
#' @param schedule Character; the gradient-step decay schedule -- how the
#'   per-iteration scale grows with the iteration index \eqn{i}. The coefficient
#'   update is \eqn{W \leftarrow W - \texttt{step\_size}\, \nabla /
#'   (\texttt{scale} \cdot \|\nabla\|)}, so a faster-growing \code{scale} means a
#'   faster-decaying effective step. One of:
#'   \describe{
#'     \item{\code{"sqrt"}}{(default) \eqn{\texttt{scale} = \sqrt{i}}, effective
#'       step \eqn{\propto 1/\sqrt{i}} -- the classic subgradient schedule and
#'       the historical hardcoded behaviour.}
#'     \item{\code{"log"}}{\eqn{\texttt{scale} = \log(i) + 1} (slow decay).}
#'     \item{\code{"const"}}{\eqn{\texttt{scale} = 1} (no decay).}
#'     \item{\code{"linear"}}{\eqn{\texttt{scale} = i}, effective step
#'       \eqn{\propto 1/i} (Robbins-Monro).}
#'   }
#'   The default \code{"sqrt"} is retained for backward compatibility: it
#'   reproduces previously published results exactly. An A/B study shipped with
#'   the package (\code{validation/ab_s02_gradient_schedule.R}) found that
#'   \code{"sqrt"} reaches the \emph{correct} optimum but converges to it roughly
#'   12--17x slower than \code{"log"} and \code{"const"} (which reach the same
#'   objective in a small fraction of the iterations), while \code{"linear"}
#'   decays too aggressively and can stall short of the optimum. Prefer
#'   \code{"log"} or \code{"const"} when convergence speed matters; note that
#'   changing the schedule changes the optimisation path and therefore the exact
#'   coefficients, so it should not be altered when reproducing published runs.
#' @param survival a logical value indicating whether survival analysis should be performed, default is FALSE
#' @param censored.mat a matrix indicating whether observations are censored, used only if survival=TRUE
#' @param tol Convergence tolerance. The algorithm stops when the
#'   improvement in the objective function is less than \code{tol}.
#'   Default: \code{1e-5}.
#' @return a list containing the following elements:
#' \item{B}{a matrix of regression coefficients with dimensions p x numTasks}
#' \item{obj}{the final objective function when algorithms stops}
#' \item{imax}{the number of iterations}
#'
#' @export
#'
#' @references Kevin Vervier, Pierre Mahé, Alexandre d’Aspremont, Jean-Baptiste Veyrieras, Jean-Philippe Vert (2014). On learning matrices with orthogonal columns or disjoint supports. https://hal.science/hal-00985654/file/learningDisjointSupports.pdf
#' @examples
#' # Regression with orthogonal columns
#' set.seed(42)
#' n <- 100; p <- 10; n_tasks <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' W_true <- qr.Q(qr(matrix(rnorm(p * n_tasks), p, n_tasks)))
#' Y <- X %*% W_true + matrix(rnorm(n * n_tasks), n) * 0.1
#' K <- matrix(1, n_tasks, n_tasks)
#' diag(K) <- 0.5
#' fit <- orthoMTL(X, Y, lambda = 1e-3, K = K, disjoint = FALSE)
#' fit$B           # coefficient matrix
#' fit$converged   # did optimisation converge?
#'
# TODO(v1.1): Loss is normalised by nrow(X) but penalty terms are not.
#   Effective regularisation strength depends on sample size.
#   Investigate impact on SOLAR-1 results and simulated vignette.

# RESOLVED (S-02, 2026-06): the gradient scaling is now the `schedule`
#   parameter (default "sqrt" reproduces the historical hardcoded behaviour).
#   The A/B in validation/ab_s02_gradient_schedule.R quantifies the trade-off:
#   "log"/"const" reach the same optimum far faster; "linear" can stall.

# TODO(v1.1): stop_no_improve default (100) may be insufficient for
#   high-dimensional or survival problems. Investigate on SOLAR-1 data.

# TODO(v1.1): W_0 initialisation differs between disjoint (random) and
#   non-disjoint (zero). This creates asymmetric seed-dependence.
#   Investigate whether consistent initialisation changes results.

orthoMTL <- function(X, Y, lambda = 1,
                     step_size = 0.1, tol = 1e-5,
                     stop_no_improve = 100, max_iter = 1e+06,
                     W_0 = NULL, seed = 42,
                     K = NULL, disjoint = FALSE, logistic = FALSE,
                     alpha = 0,
                     schedule = c("sqrt", "log", "const", "linear"),
                     survival = FALSE, censored.mat = NULL,
                     verbose = 0){
  # gradient-step decay schedule (S-02); default "sqrt" = historical behaviour
  schedule <- match.arg(schedule)
  # initiate the random seed
  set.seed(seed)

  # extract problem dimensions
  if(is.null(X) | is.null(Y)) stop('Training data X and Y need to be provided \n')
  m = nrow(X)
  p = ncol(X)
  numTasks = ncol(Y)

  # regularization parameters test
  if(lambda <0) stop('Regularization parameter lambda needs to be positive value \n')

  # if the K matrix is not defined, set it as identity matrix
  if (is.null(K)) K = diag(1, numTasks)
  # if no initial W_0 was provided, generate a random one.
  if(is.null(W_0)){
    # in the non-disjoint case, we observed that results were stable when starting with the zero matrix
    if (!disjoint) {
      W_k = matrix(0, nrow = p, ncol = numTasks)
    }
    else {
      # in the disjoint case, we initiate with a random matrix with positive values
      W_k = matrix(rnorm(p * numTasks), ncol = numTasks)
      W_k = abs(W_k)
    }
  }else{
    W_k = W_0
  }
  # store current model
  W = W_k
  # in the disjoint case, a complementary matrix V is created
  if (disjoint) {
    V_k = W_k
    V <- V_k
  }

  # elastic-net mixing parameter must lie in [0, 1]
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("Mixing parameter alpha must be a single number in [0, 1]", call. = FALSE)
  }

  # initiate variables tracking progress
  new <- Inf # new objective value
  no_improv = 0 # number of steps without improvement
  i = 0 # step counter

  # print header
  if (verbose > 0) cat("Step \t ObjFun\t NonImproving\n", sep = "")

  # main loop
  while ((no_improv < stop_no_improve) && (i < max_iter)) {
    i <- i + 1
    # print update every 1000 steps
    if (verbose > 0 & i%%1000 == 0) {
      cat(i, "\t", new, "\t", no_improv, "\n", sep = "")
      idx = which(W_k != 0, arr.ind = TRUE)
      cat("Current sparse support:", length(W_k) - nrow(idx),
          "\n")
    }

    ####################
    # loss calculation #
    #------------------#

    # loss pre-computation
    # precompute the linear predictor once and reuse in loss and gradient
    XW_k = X %*% W_k
    if (!logistic) {
      LS = XW_k - Y #MSE
    }
    else {
      # Numerically stable softplus: log(1+exp(-z)) = max(-z,0) + log1p(exp(-|z|))
      # where z = Y * (X W). Avoids exp() overflow for large |z|. (REG-02)
      score = Y * XW_k
      LS = pmax(-score, 0) + log1p(exp(-abs(score))) # Logistic
    }

    # for survival analysis, the loss calculation is restricted to non-censored data
    if(survival){
      if(!is.null(censored.mat)){
        LS[censored.mat == 0] = 0
      }else{stop("A matrix with censoring information needs to be provided\n")}
    }

    #compute a scale for gradient descent (S-02: `schedule` controls the decay)
    scale = switch(schedule,
                   sqrt   = sqrt(i),
                   log    = log(i) + 1,
                   const  = 1,
                   linear = i)

    # orthogonality penalty matrix (V^T V for disjoint, W^T W otherwise)
    if (disjoint) {
      PEN = crossprod(V_k)
      W_or_V = V_k
    }else {
      PEN = crossprod(W_k)
      W_or_V = W_k
    }

    # elastic-net penalty:
    #   lambda * [ (1 - alpha)/2 * Omega_K(W)^2  +  alpha * ||W||_1 ]
    # alpha = 0 -> pure orthogonality penalty; alpha = 1 -> pure Lasso.
    pen_obj = (1 - alpha) / 2 * lambda * sum(abs(PEN) * K) +
      alpha * lambda * sum(abs(W_or_V))

    # total objective calculation
    if (!logistic) {
      tmp = 0.5 * sum(LS^2) / nrow(X) + pen_obj
    }else {
      tmp = sum(LS) / nrow(X) + pen_obj
    }

    # Previously: abs(tmp - new) > 10^-5
    if (tmp < new && abs(tmp - new) > tol) {
      # if improvment, update all terms
      no_improv <- 0
      new <- tmp
      W <- W_k
      if (disjoint)
        V <- V_k
    }
    else { #if no improvement, add one
      no_improv <- no_improv + 1
    }

    ########################
    # gradient computation #
    #----------------------#

    if (disjoint) {
      # gradient of the data-fit term w.r.t. W
      if (!logistic) {
        gradientW <- crossprod(X, LS)/nrow(X)
      }else {
        # Correct subgradient for y in {-1,+1}: -X^T(y/(1+exp(y*Xw)))/m (REG-01)
        gradientW <- -crossprod(X, Y/(1 + exp(pmin(Y * XW_k, 500))))/nrow(X)
      }
      # subgradient of the elastic-net penalty w.r.t. V
      gradientV <- (1 - alpha) * lambda * V_k %*% (sign(PEN) * K) +
        alpha * lambda * sign(V_k)
      norm_gradient <- sqrt(sum((gradientW + gradientV)^2))
    }else {
      # subgradient of data-fit + elastic-net penalty w.r.t. W
      if (!logistic) {
        gradientW <- crossprod(X, LS)/nrow(X) +
          (1 - alpha) * lambda * W_k %*% (sign(PEN) * K) +
          alpha * lambda * sign(W_k)
      }else {
        gradientW <- -crossprod(X, Y/(1 + exp(pmin(Y * XW_k, 500))))/nrow(X) +
          (1 - alpha) * lambda * W_k %*% (sign(PEN) * K) +
          alpha * lambda * sign(W_k)
      }
      norm_gradient <- sqrt(sum((gradientW)^2))
    }

    ############################################
    # update coefficients via gradient descent #
    #------------------------------------------#

    if (norm_gradient > 0) {
      scale_grad = scale * norm_gradient

      W_k <- W_k - step_size * gradientW/scale_grad
      if (disjoint)
        V_k <- V_k - step_size * gradientV/scale_grad

      if (disjoint) {
        projection <- proj_disjoint(w = W_k, v = V_k)
        W_k <- projection$w
        V_k <- projection$v
      }

    }else { #if norm is zero, we stop
      no_improv <- stop_no_improve + 1
    }
  }

  # create final output
  imax = i
  if(imax == max_iter) warning('Optimization reached the maximal number of iterations without converging \n
                               Please consider increasing the number of iterations or the step size \n')
  best_obj <- new
  w_star <- W
  # Build S3 return object
  structure(
    list(
      B              = w_star,
      W              = w_star,  # orthopen backward compatibility alias
      obj            = new,
      imax           = imax,
      converged      = imax < max_iter,
      call           = match.call(),
      hyperparameters = list(
        lambda          = lambda,
        alpha           = alpha,
        schedule        = schedule,
        step_size       = step_size,
        tol             = tol,
        stop_no_improve = stop_no_improve,
        max_iter        = max_iter,
        seed            = seed,
        disjoint        = disjoint,
        logistic        = logistic,
        survival        = survival
      ),
      K              = K,
      n_tasks        = numTasks,
      n_features     = p,
      n_obs          = m,
      feature_names  = colnames(X),
      task_names     = colnames(Y)
    ),
    class = "orthoMTL"
  )
}
