#' Bootstrap Inference for orthoMTL Coefficients
#'
#' Estimates coefficient variability and statistical relevance by
#' comparing bootstrapped models (real signal) against null models
#' (permuted outcomes). This two-pronged approach answers: (1) how
#' stable is each coefficient across resamples, and (2) is each
#' coefficient distinguishable from what would be obtained by chance.
#'
#' @param X A numeric matrix of predictor variables with dimensions
#'   \code{n x p}.
#' @param Y A numeric matrix of response labels with dimensions
#'   \code{n x numTasks}. May contain \code{NA} for censored
#'   observations.
#' @param lambda Regularisation parameter for the orthogonal penalty.
#' @param alpha Elastic-net mixing parameter in \eqn{[0, 1]}, passed
#'   through to \code{\link{orthoMTL}}. Default: \code{0}.
#' @param step_size Step size for gradient descent. Default: \code{0.1}.
#' @param K Constraint matrix of dimensions
#'   \code{numTasks x numTasks}. Default: \code{NULL} (identity).
#' @param disjoint Logical. Enforce disjoint supports? Default:
#'   \code{FALSE}.
#' @param survival Logical. Use censored survival loss? Default:
#'   \code{FALSE}.
#' @param censored.mat A numeric indicator matrix of dimensions
#'   \code{n x numTasks}. Required when \code{survival = TRUE}.
#' @param n_repeats Number of bootstrap/permutation repeats. Default:
#'   \code{100}.
#' @param n_cores Number of cores for parallel execution. Default:
#'   \code{2}.
#' @param seed Base random seed. Each repeat uses \code{set.seed(i)}
#'   for \code{i} in \code{1:n_repeats} to ensure reproducibility.
#'   Default: \code{42} (unused directly; individual repeats use
#'   their index as seed).
#' @param verbose Logical. Print progress information? Default:
#'   \code{TRUE}.
#'
#' @return An object of class \code{"bootstrap_orthoMTL"} containing:
#'   \describe{
#'     \item{results}{A tidy \code{data.frame} with columns \code{id}
#'       (feature name), \code{time} (task/threshold), \code{coeff}
#'       (coefficient value), \code{group} (\code{"real"} or
#'       \code{"null"}), and \code{repeat_id} (integer).}
#'     \item{coefficients_real}{A list of \code{n_repeats} coefficient
#'       matrices (each \code{p x numTasks}).}
#'     \item{coefficients_null}{A list of \code{n_repeats} coefficient
#'       matrices from permuted outcomes.}
#'     \item{obj_real}{Numeric vector of final objective values for
#'       real bootstrap models.}
#'     \item{obj_null}{Numeric vector of final objective values for
#'       null permutation models.}
#'     \item{n_repeats}{Number of repeats.}
#'     \item{n_features}{Number of features.}
#'     \item{n_tasks}{Number of tasks.}
#'     \item{feature_names}{Character vector of feature names.}
#'     \item{task_names}{Character vector of task names.}
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' \strong{Real bootstrap}: For each repeat, rows of \code{X}, \code{Y},
#' and \code{censored.mat} are resampled \emph{with replacement}. The
#' model is refit with identical hyperparameters. This produces a
#' distribution of coefficient values reflecting estimation variability.
#'
#' \strong{Null permutation}: For each repeat, rows of \code{Y} and
#' \code{censored.mat} are permuted \emph{without replacement} while
#' \code{X} remains fixed. This breaks the association between features
#' and outcomes, producing a null distribution of coefficients.
#'
#' Comparing real vs null distributions for each feature and task
#' indicates whether observed coefficients are distinguishable from
#' noise. Visualise with \code{\link{plot_bootstrap}}.
#'
#' @seealso \code{\link{orthoMTL}}, \code{\link{plot_bootstrap}}
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 30; p <- 5; n_tasks <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' SurvTime <- rexp(n, rate = 0.1)
#' Event <- rbinom(n, 1, 0.7)
#' thresholds <- c(4, 6, 10)
#'
#' Y <- create_longitudinal_labels(SurvTime, Event, thresholds)
#' W <- create_indicator_matrix(Y)
#' K <- create_constraint_matrix(n_tasks)
#'
#' boot_res <- bootstrap_orthoMTL(
#'   X = X, Y = Y, lambda = 1e-3, step_size = 0.1,
#'   K = K, survival = TRUE, censored.mat = W,
#'   n_repeats = 5, n_cores = 1, verbose = FALSE
#' )
#'
#' print(boot_res)
#' head(boot_res$results)
#' }
bootstrap_orthoMTL <- function(X, Y,
                               lambda = 1,
                               alpha = 0,
                               step_size = 0.1,
                               K = NULL,
                               disjoint = FALSE,
                               survival = FALSE,
                               censored.mat = NULL,
                               n_repeats = 100,
                               n_cores = 2,
                               seed = 42,
                               verbose = TRUE) {

  cl <- match.call()

  # ---------------------------
  # Input validation
  # ---------------------------
  if (is.null(X) || is.null(Y)) {
    stop("X and Y must be provided.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows.", call. = FALSE)
  }
  if (survival && is.null(censored.mat)) {
    stop("censored.mat must be provided when survival = TRUE.",
         call. = FALSE)
  }
  if (n_repeats < 1) {
    stop("n_repeats must be at least 1.", call. = FALSE)
  }

  n <- nrow(X)
  p <- ncol(X)
  numTasks <- ncol(Y)
  feature_names <- colnames(X)
  task_names <- colnames(Y)

  if (is.null(feature_names)) {
    feature_names <- paste0("V", seq_len(p))
  }
  if (is.null(task_names)) {
    task_names <- paste0("T", seq_len(numTasks))
  }

  # ---------------------------
  # Parallel setup
  # ---------------------------
  if (n_cores > 1L) {
    cl_parallel <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl_parallel)
    # Workers are fresh R sessions (PSOCK on Windows); ensure they can locate
    # orthoMTL even when it lives in a temporary library (e.g. during R CMD
    # check or load_all) by inheriting the parent's library paths.
    parallel::clusterCall(cl_parallel, function(p) .libPaths(p), .libPaths())
    on.exit(parallel::stopCluster(cl_parallel), add = TRUE)
  } else {
    # Single core: run sequentially in-process. Avoids spawning a PSOCK worker
    # that cannot library(orthoMTL) under R CMD check (the Windows-only
    # vignette build failure).
    foreach::registerDoSEQ()
  }

  # ---------------------------
  # Phase 1: Real bootstrap
  # ---------------------------
  if (verbose) {
    cat("Bootstrap inference for orthoMTL\n")
    cat("  Features:", p, "| Tasks:", numTasks, "| Repeats:", n_repeats, "\n")
    cat("  Cores:", n_cores, "\n")
    cat("  Phase 1/2: Real bootstrap (resampling with replacement)...\n")
  }

  real_results <- foreach::foreach(
    i = seq_len(n_repeats),
    .packages = "orthoMTL"
  ) %dopar% {

    set.seed(i)
    idx <- sample(x = seq_len(n), size = n, replace = TRUE)

    fit <- orthoMTL(
      X            = X[idx, , drop = FALSE],
      Y            = Y[idx, , drop = FALSE],
      lambda       = lambda,
      alpha        = alpha,
      step_size    = step_size,
      K            = K,
      disjoint     = disjoint,
      survival     = survival,
      censored.mat = if (survival) censored.mat[idx, , drop = FALSE] else NULL,
      seed         = i,
      verbose      = 0
    )

    Bt <- fit$B
    colnames(Bt) <- task_names
    rownames(Bt) <- feature_names

    list(B = Bt, obj = fit$obj)
  }

  # Unpack real results
  coefficients_real <- lapply(real_results, `[[`, "B")
  obj_real <- vapply(real_results, `[[`, numeric(1), "obj")

  if (verbose) cat("  Done.\n")

  # ---------------------------
  # Phase 2: Null permutation
  # ---------------------------
  if (verbose) {
    cat("  Phase 2/2: Null permutation (shuffled outcomes)...\n")
  }

  null_results <- foreach::foreach(
    i = seq_len(n_repeats),
    .packages = "orthoMTL"
  ) %dopar% {

    set.seed(i)
    idx <- sample(x = seq_len(n), size = n, replace = FALSE)

    fit <- orthoMTL(
      X            = X,
      Y            = Y[idx, , drop = FALSE],
      lambda       = lambda,
      alpha        = alpha,
      step_size    = step_size,
      K            = K,
      disjoint     = disjoint,
      survival     = survival,
      censored.mat = if (survival) censored.mat[idx, , drop = FALSE] else NULL,
      seed         = i,
      verbose      = 0
    )

    Bt <- fit$B
    colnames(Bt) <- task_names
    rownames(Bt) <- feature_names

    list(B = Bt, obj = fit$obj)
  }

  # Unpack null results
  coefficients_null <- lapply(null_results, `[[`, "B")
  obj_null <- vapply(null_results, `[[`, numeric(1), "obj")

  if (verbose) cat("  Done.\n")

  # ---------------------------
  # Build tidy results data.frame
  # ---------------------------
  if (verbose) cat("  Assembling results...\n")

  results_list <- vector("list", length(feature_names))

  for (f_idx in seq_along(feature_names)) {
    feature <- feature_names[f_idx]

    # Real coefficients: n_repeats x numTasks matrix
    real_mat <- do.call("rbind", lapply(coefficients_real,
                                        function(B) B[feature, ]))

    real_df <- data.frame(
      id        = feature,
      time      = rep(task_names, each = n_repeats),
      coeff     = as.vector(real_mat),
      group     = "real",
      repeat_id = rep(seq_len(n_repeats), times = numTasks),
      stringsAsFactors = FALSE
    )

    # Null coefficients: n_repeats x numTasks matrix
    null_mat <- do.call("rbind", lapply(coefficients_null,
                                        function(B) B[feature, ]))

    null_df <- data.frame(
      id        = feature,
      time      = rep(task_names, each = n_repeats),
      coeff     = as.vector(null_mat),
      group     = "null",
      repeat_id = rep(seq_len(n_repeats), times = numTasks),
      stringsAsFactors = FALSE
    )

    results_list[[f_idx]] <- rbind(real_df, null_df)
  }

  results <- do.call("rbind", results_list)
  rownames(results) <- NULL

  if (verbose) cat("  Complete.\n")

  # ---------------------------
  # Return S3 object
  # ---------------------------
  structure(
    list(
      results           = results,
      coefficients_real = coefficients_real,
      coefficients_null = coefficients_null,
      obj_real          = obj_real,
      obj_null          = obj_null,
      n_repeats         = n_repeats,
      n_features        = p,
      n_tasks           = numTasks,
      feature_names     = feature_names,
      task_names        = task_names,
      call              = cl
    ),
    class = "bootstrap_orthoMTL"
  )
}
