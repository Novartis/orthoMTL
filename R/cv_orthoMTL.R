#' Cross-Validation for orthoMTL Hyperparameter Selection
#'
#' Performs a parallelised grid search over hyperparameters for
#' \code{\link{orthoMTL}}, evaluating each configuration via
#' cross-validated concordance index. Returns the best configuration
#' without retraining a final model (that is the caller's responsibility).
#'
#' @param X.train A numeric matrix of training features with dimensions
#'   \code{n x p}.
#' @param Y.train A numeric matrix of training labels with dimensions
#'   \code{n x numTasks}. May contain \code{NA} for censored observations.
#' @param W.train A numeric indicator matrix of dimensions
#'   \code{n x numTasks}, where \code{1} = observed and \code{0} = censored.
#'   Required when \code{survival = TRUE}.
#' @param K A square constraint matrix of dimensions
#'   \code{numTasks x numTasks}. The diagonal will be overridden by
#'   values in \code{diag_vals} during the search.
#' @param lambdas A numeric vector of regularisation parameters to search.
#' @param alphas A numeric vector of elastic-net mixing parameters in
#'   \eqn{[0, 1]} to search. Default: \code{0} (pure orthogonality
#'   penalty, no sparsity).
#' @param stepsizes A numeric vector of gradient descent step sizes
#'   to search.
#' @param diag_vals A numeric vector of diagonal values for the
#'   constraint matrix \code{K} to search.
#' @param survival Logical. Use censored survival loss? Default:
#'   \code{TRUE}.
#' @param logistic Logical. Fit logistic (classification) models?
#'   Passed through to \code{\link{orthoMTL}}. Default: \code{FALSE}.
#' @param metric Character; the scoring metric maximised/minimised over
#'   the grid, or \code{NULL} (default) to choose automatically by mode:
#'   \code{"cindex"} when \code{survival = TRUE}, \code{"auc"} when
#'   \code{logistic = TRUE}, otherwise \code{"rmse"}. Supported values:
#'   \code{"cindex"}, \code{"auc"}, \code{"accuracy"} (higher is better),
#'   \code{"rmse"} (lower is better), \code{"r2"} (higher is better).
#' @param disjoint Logical. Enforce disjoint supports? Default:
#'   \code{FALSE}.
#' @param schedule Character; the gradient-step decay schedule passed to
#'   \code{\link{orthoMTL}} for every fit. One of \code{"sqrt"} (default),
#'   \code{"log"}, \code{"const"}, \code{"linear"}. \code{"log"} or
#'   \code{"const"} typically reach the same optimum in far fewer iterations
#'   than the default \code{"sqrt"}, which can noticeably speed up the grid
#'   search; see the \code{schedule} argument of \code{\link{orthoMTL}} for the
#'   trade-offs. Applied uniformly to all configurations (it is not part of the
#'   tuning grid).
#' @param folds An integer vector of length \code{n} assigning each
#'   training observation to a fold. If \code{NULL}, a 5-fold
#'   assignment is generated with a warning.
#' @param n_cores Integer. Number of cores for parallel execution.
#'   Default: \code{2}.
#' @param seed Optional integer random seed for reproducibility. Default:
#'   \code{NULL} (no seed is set; the ambient RNG state is used as-is).
#'   Used both for auto-generated fold assignment and for each
#'   \code{\link{orthoMTL}} fit in the grid.
#' @param verbose Logical. Print progress information? Default:
#'   \code{TRUE}.
#'
#' @return An object of class \code{"cv_orthoMTL"} containing:
#'   \describe{
#'     \item{best}{A list with the best hyperparameters: \code{lambda},
#'       \code{alpha}, \code{stepsize}, \code{diag_val}, and the
#'       corresponding \code{cv_score}.}
#'     \item{results}{A \code{data.frame} of all configurations with
#'       their mean CV C-index, sorted descending by \code{cv_score}.}
#'     \item{folds}{The fold assignment vector used.}
#'     \item{n_configs}{Total number of configurations tested.}
#'     \item{n_folds}{Number of unique folds.}
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' The grid is constructed as the full Cartesian product of
#' \code{lambdas}, \code{alphas}, \code{stepsizes}, and
#' \code{diag_vals}. Each configuration is evaluated independently
#' in parallel across cores. Within each configuration, folds are
#' evaluated sequentially and the per-fold C-indices are averaged.
#'
#' The best configuration is selected by joint maximisation of the
#' mean CV C-index over the entire flattened grid (not greedy
#' sequential search).
#'
#' This function does \emph{not} retrain a final model. Use the
#' returned hyperparameters to train via \code{\link{orthoMTL}}.
#'
#' @seealso \code{\link{orthoMTL}}, \code{\link{cindex_mtl}}
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
#' n <- 50; p <- 5; n_tasks <- 3
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
#' folds <- rep(1:2, length.out = n)
#'
#' cv_res <- cv_orthoMTL(
#'   X.train = X, Y.train = Y, W.train = W, K = K,
#'   lambdas = c(1e-3, 1e-2), alphas = 0,
#'   stepsizes = c(0.1), diag_vals = c(0.5, 1),
#'   survival = TRUE, disjoint = FALSE,
#'   folds = folds, n_cores = 1, seed = 42, verbose = FALSE
#' )
#'
#' print(cv_res)
#' cv_res$best
#' }
cv_orthoMTL <- function(X.train, Y.train, W.train = NULL,
                        K = NULL,
                        lambdas = c(1e-3, 1e-2),
                        alphas = 0,
                        stepsizes = c(0.1, 0.5),
                        diag_vals = c(0.5, 1),
                        survival = TRUE,
                        logistic = FALSE,
                        metric = NULL,
                        disjoint = FALSE,
                        schedule = c("sqrt", "log", "const", "linear"),
                        folds = NULL,
                        n_cores = 2,
                        seed = NULL,
                        verbose = TRUE) {

  cl <- match.call()
  schedule <- match.arg(schedule)

  # ---------------------------
  # Input validation
  # ---------------------------
  if (is.null(X.train) || is.null(Y.train)) {
    stop("X.train and Y.train must be provided.", call. = FALSE)
  }
  if (nrow(X.train) != nrow(Y.train)) {
    stop("X.train and Y.train must have the same number of rows.",
         call. = FALSE)
  }
  if (survival && is.null(W.train)) {
    stop("W.train (censoring indicator) must be provided when ",
         "survival = TRUE.", call. = FALSE)
  }
  if (survival && !is.null(W.train)) {
    if (nrow(W.train) != nrow(X.train) || ncol(W.train) != ncol(Y.train)) {
      stop("W.train dimensions must match Y.train.", call. = FALSE)
    }
  }

  numTasks <- ncol(Y.train)

  # Default K

  if (is.null(K)) {
    K <- diag(1, numTasks)
  }

  # ---------------------------
  # Fold assignment
  # ---------------------------
  if (is.null(folds)) {
    warning("Cross-validation fold assignment not provided. ",
            "Generating 5-fold assignment.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)
    folds <- sample(1:5, nrow(X.train), replace = TRUE)
  }
  if (length(folds) != nrow(X.train)) {
    stop("Length of 'folds' (", length(folds), ") must equal nrow(X.train) (",
         nrow(X.train), ").", call. = FALSE)
  }
  unique_folds <- sort(unique(folds))
  n_folds <- length(unique_folds)

  # ---------------------------
  # Elastic-net mixing validation
  # ---------------------------
  if (any(alphas < 0 | alphas > 1)) {
    stop("All 'alphas' must lie in [0, 1].", call. = FALSE)
  }

  # ---------------------------
  # Resolve scoring metric (mode-aware default)
  # ---------------------------
  scorer <- .resolve_cv_metric(metric, survival = survival,
                               logistic = logistic)

  # ---------------------------
  # Build configuration grid
  # ---------------------------
  config_grid <- expand.grid(
    lambda   = lambdas,
    alpha    = alphas,
    stepsize = stepsizes,
    diag_val = diag_vals,
    stringsAsFactors = FALSE
  )
  n_configs <- nrow(config_grid)

  if (verbose) {
    total_fits <- n_configs * n_folds
    cat("Cross-validation for orthoMTL\n")
    cat("  Configurations:", n_configs, "\n")
    cat("  Folds:", n_folds, "\n")
    cat("  Total model fits:", total_fits, "\n")
    cat("  Cores:", n_cores, "\n")
    cat("  Running...\n")
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
    # Ensure cluster is stopped on exit (even on error)
    on.exit(parallel::stopCluster(cl_parallel), add = TRUE)
  } else {
    # Single core: run sequentially in-process. Avoids spawning a PSOCK worker
    # that cannot library(orthoMTL) under R CMD check (the Windows-only
    # vignette build failure).
    foreach::registerDoSEQ()
  }

  # ---------------------------
  # Main CV loop (parallel over configs)
  # ---------------------------
  cv_scores <- foreach::foreach(
    cfg_idx = seq_len(n_configs),
    .combine  = "c",
    .packages = "orthoMTL"
  ) %dopar% {

    cfg <- config_grid[cfg_idx, ]

    # Local copy of K with this config's diagonal
    K_local <- K
    diag(K_local) <- cfg$diag_val

    # Sequential fold loop
    fold_scores <- numeric(n_folds)

    for (f_idx in seq_along(unique_folds)) {
      fold_id <- unique_folds[f_idx]
      idx_train <- which(folds != fold_id)
      idx_val   <- which(folds == fold_id)

      # Fit model on training folds
      fit <- orthoMTL(
        X            = X.train[idx_train, , drop = FALSE],
        Y            = Y.train[idx_train, , drop = FALSE],
        lambda       = cfg$lambda,
        alpha        = cfg$alpha,
        step_size    = cfg$stepsize,
        K            = K_local,
        disjoint     = disjoint,
        survival     = survival,
        logistic     = logistic,
        schedule     = schedule,
        censored.mat = if (survival) W.train[idx_train, , drop = FALSE] else NULL,
        seed         = seed,
        verbose      = 0
      )

      # Predict on validation fold
      X_val <- X.train[idx_val, , drop = FALSE]
      Y_val <- Y.train[idx_val, , drop = FALSE]

      pred_val <- X_val %*% fit$B

      # Apply monotonicity projection for survival
      if (survival) {
        pred_val <- t(apply(pred_val, 1, nnmaxheap_C))
      }

      # Score on validation fold with the resolved metric
      fold_scores[f_idx] <- scorer$fn(Y_val, pred_val)
    }

    # Return mean CV score for this config
    mean(fold_scores, na.rm = TRUE)
  }

  # ---------------------------
  # Collect results
  # ---------------------------
  config_grid$cv_score <- cv_scores

  # Joint optimum over flattened grid (direction depends on the metric)
  best_idx <- if (scorer$maximize) {
    which.max(config_grid$cv_score)
  } else {
    which.min(config_grid$cv_score)
  }
  best_cfg <- config_grid[best_idx, ]

  if (verbose) {
    cat("  Done.\n")
    cat("  Best CV ", scorer$name, ": ",
        format(best_cfg$cv_score, digits = 4), "\n", sep = "")
  }

  # Check for failed fits. NA always indicates a failure; an exact 0 is
  # suspicious only for the "higher-is-better" metrics (a 0 RMSE is perfect).
  n_failed <- sum(is.na(cv_scores) |
                    (scorer$maximize & cv_scores == 0))
  if (n_failed > 0) {
    warning(n_failed, " of ", n_configs,
            " configurations returned NA or a degenerate score. ",
            "This may indicate failed convergence or degenerate folds.",
            call. = FALSE)
  }

  # Sort results best-first (direction depends on the metric)
  config_grid <- config_grid[order(config_grid$cv_score,
                                   decreasing = scorer$maximize), ]
  rownames(config_grid) <- NULL

  # ---------------------------
  # Build return object
  # ---------------------------
  structure(
    list(
      best = list(
        lambda   = best_cfg$lambda,
        alpha    = best_cfg$alpha,
        stepsize = best_cfg$stepsize,
        diag_val = best_cfg$diag_val,
        cv_score = best_cfg$cv_score
      ),
      results   = config_grid,
      folds     = folds,
      n_configs = n_configs,
      n_folds   = n_folds,
      metric    = scorer$name,
      call      = cl
    ),
    class = "cv_orthoMTL"
  )
}


# ---------------------------------------------------------------------
# Internal: resolve the CV scoring metric (GEN-02)
#
# Returns a list with the metric name, the scoring function fn(true,
# pred), and `maximize` (TRUE if higher is better). When `metric` is
# NULL the default is chosen from the fitting mode.
# ---------------------------------------------------------------------
.resolve_cv_metric <- function(metric, survival, logistic) {
  if (is.null(metric)) {
    metric <- if (survival) {
      "cindex"
    } else if (logistic) {
      "auc"
    } else {
      "rmse"
    }
  }
  metric <- match.arg(metric,
                      c("cindex", "auc", "accuracy", "rmse", "r2"))

  switch(metric,
    cindex   = list(name = "cindex",   fn = cindex_mtl,   maximize = TRUE),
    auc      = list(name = "auc",      fn = auc_mtl,      maximize = TRUE),
    accuracy = list(name = "accuracy", fn = accuracy_mtl, maximize = TRUE),
    r2       = list(name = "r2",       fn = r2_mtl,       maximize = TRUE),
    rmse     = list(name = "rmse",     fn = rmse_mtl,     maximize = FALSE)
  )
}
