#' Simulate Multi-Task Data with Time-Varying Effects
#'
#' Generates a realistic simulated dataset with binary and continuous
#' features and a known ground-truth coefficient structure. The
#' \code{mode} argument selects the response type: piecewise-exponential
#' \strong{survival} times (default), multi-task \strong{regression}
#' targets, or multi-task binary \strong{classification} labels. Designed
#' for demonstrating and testing \code{\link{orthoMTL}}.
#'
#' @param mode Character; the response type to generate. One of
#'   \code{"survival"} (default), \code{"regression"}, or
#'   \code{"classification"}. The feature matrix and ground-truth
#'   coefficients are generated identically across modes; only the
#'   response differs:
#'   \describe{
#'     \item{survival}{piecewise-exponential \code{SurvTime}/\code{Event}
#'       (the historical behaviour).}
#'     \item{regression}{\code{Y = X \%*\% beta + N(0, noise_sd^2)},
#'       returned as the \code{n x length(thresholds)} matrix \code{Y}.}
#'     \item{classification}{binary labels in \eqn{\{-1, +1\}} sampled
#'       from \eqn{P(Y = +1) = 1 / (1 + e^{-X beta})}, returned as
#'       \code{Y}. This matches the label encoding \code{orthoMTL}
#'       optimises with \code{logistic = TRUE}.}
#'   }
#' @param noise_sd Standard deviation of the Gaussian noise added in
#'   \code{mode = "regression"}. Ignored otherwise. Default: \code{1}.
#' @param n Number of patients. Default: \code{200}.
#' @param p Number of features excluding the treatment column.
#'   Default: \code{30}.
#' @param n_signals Number of features with true non-zero effects.
#'   Default: \code{5}. Must be \code{<= p}.
#' @param n_continuous Number of continuous features (placed first
#'   in the feature matrix). Default: \code{1}. Remaining features
#'   are binary.
#' @param thresholds Numeric vector of time thresholds defining the
#'   task structure. Default: \code{c(4, 6, 10, 15)}.
#' @param censoring_max Maximum censoring time. Censoring times are
#'   drawn from \code{Unif(0, censoring_max)}. Default: \code{25}.
#' @param baseline_hazard Baseline hazard rate per interval. Default:
#'   \code{0.05}.
#' @param effect_strength Multiplier controlling the magnitude of
#'   feature effects. Default: \code{0.5}.
#' @param treatment_effect Effect of treatment on the log-hazard
#'   (negative = protective). Default: \code{-0.3}.
#' @param seed Optional random seed for reproducibility. Default:
#'   \code{NULL} (no seed is set; the ambient RNG state is used as-is).
#'   Pass an integer to make the simulated data reproducible.
#'
#' @return An object of class \code{"simulated_mtl"} containing:
#'   \describe{
#'     \item{mode}{The response type generated.}
#'     \item{Y}{For \code{mode = "regression"} / \code{"classification"},
#'       the \code{n x length(thresholds)} response matrix. \code{NULL}
#'       in survival mode (use \code{SurvTime}/\code{Event} instead).}
#'     \item{SurvTime, Event}{Survival mode only (\code{NULL} otherwise):
#'       observed times and event indicators.}
#'     \item{X}{Numeric matrix of dimensions \code{n x (p + 1)}.
#'       Columns include \code{n_continuous} continuous features
#'       (standard normal), \code{p - n_continuous} binary features
#'       (prevalences drawn from \code{Beta(2, 10)}), and a
#'       \code{treatment} column (balanced 1:1).}
#'     \item{SurvTime}{Numeric vector of observed survival times.}
#'     \item{Event}{Binary vector: \code{1} = event observed,
#'       \code{0} = censored.}
#'     \item{treatment}{Binary vector (also present as last column
#'       of \code{X}).}
#'     \item{feature_names}{Character vector of column names of
#'       \code{X}.}
#'     \item{ground_truth}{A list containing:
#'       \describe{
#'         \item{coefficients}{Matrix of true coefficients with
#'           dimensions \code{(p + 1) x length(thresholds)}.}
#'         \item{signal_features}{Names of features with non-zero
#'           effects.}
#'         \item{null_features}{Names of features with zero effects.}
#'         \item{effect_types}{Named character vector mapping signal
#'           features to their temporal effect template.}
#'         \item{thresholds}{The thresholds used.}
#'         \item{baseline_hazard}{The baseline hazard used.}
#'       }
#'     }
#'     \item{n, p, n_signals, n_continuous, thresholds, seed}{Input
#'       parameters stored for reference.}
#'     \item{call}{The matched function call.}
#'   }
#'
#' @details
#' \strong{Feature structure:}
#' \itemize{
#'   \item Continuous features are drawn from \code{N(0, 1)}.
#'   \item Binary features have prevalences drawn from \code{Beta(2, 10)},
#'     producing a realistic range (~5-30\%).
#'   \item Treatment is balanced 1:1 via random assignment.
#' }
#'
#' \strong{Effect templates:}
#' Each signal feature is assigned one of five temporal patterns:
#' \itemize{
#'   \item \code{"early"}: strong effect at early thresholds, fading to
#'     zero at late.
#'   \item \code{"late"}: zero at early thresholds, emerging at late.
#'   \item \code{"constant"}: equal effect across all thresholds
#'     (detectable by standard Cox models).
#'   \item \code{"increasing"}: effect grows over time.
#'   \item \code{"decreasing"}: effect shrinks over time.
#' }
#' Templates are assigned cyclically across signal features with
#' random sign (risk-increasing or protective).
#'
#' \strong{Survival time generation:}
#' Uses a piecewise-exponential model where the hazard in each
#' interval is
#' \code{h_k(i) = baseline_hazard * exp(X[i, ] \%*\% beta[, k])}.
#' Patients progress through intervals sequentially; an event
#' occurs when the simulated time within an interval is shorter
#' than the interval width.
#'
#' \strong{Censoring:}
#' Independent of event times. \code{C ~ Unif(0, censoring_max)}.
#' Observed time = \code{min(T, C)}, event indicator = \code{T <= C}.
#'
#' @seealso \code{\link{create_longitudinal_labels}},
#'   \code{\link{orthoMTL}}
#'
#' @export
#'
#' @examples
#' # Generate simulated data
#' sim <- simulate_mtl(n = 100, p = 20, n_signals = 4, seed = 42)
#' sim
#'
#' # Inspect ground truth
#' sim$ground_truth$signal_features
#' sim$ground_truth$effect_types
#' sim$ground_truth$coefficients[sim$ground_truth$signal_features, ]
#'
#' # Use in orthoMTL workflow
#' \donttest{
#' thresholds <- c(4, 6, 10, 15)
#' Y <- create_longitudinal_labels(sim$SurvTime, sim$Event, thresholds)
#' W <- create_indicator_matrix(Y)
#' K <- create_constraint_matrix(length(thresholds))
#'
#' fit <- orthoMTL(sim$X, Y, lambda = 1e-3, K = K,
#'                 survival = TRUE, censored.mat = W)
#' summary(fit)
#' }
simulate_mtl <- function(n = 200,
                         p = 30,
                         n_signals = 5,
                         n_continuous = 1,
                         thresholds = c(4, 6, 10, 15),
                         mode = c("survival", "regression",
                                  "classification"),
                         noise_sd = 1,
                         censoring_max = 25,
                         baseline_hazard = 0.05,
                         effect_strength = 0.8,  # was 0.5
                         treatment_effect = -0.3,
                         seed = NULL) {

  cl <- match.call()
  mode <- match.arg(mode)

  # ---------------------------
  # Input validation
  # ---------------------------
  if (n < 10) stop("n must be at least 10.", call. = FALSE)
  if (p < 1) stop("p must be at least 1.", call. = FALSE)
  if (n_signals > p) {
    stop("n_signals (", n_signals, ") cannot exceed p (", p, ").",
         call. = FALSE)
  }
  if (n_continuous > p) {
    stop("n_continuous (", n_continuous, ") cannot exceed p (", p, ").",
         call. = FALSE)
  }
  if (n_signals > (p - n_continuous) + n_continuous) {
    # This is always true, but guard against n_signals > p
  }
  if (length(thresholds) < 2) {
    stop("At least 2 thresholds are required.", call. = FALSE)
  }
  if (any(diff(thresholds) <= 0)) {
    stop("Thresholds must be strictly increasing.", call. = FALSE)
  }
  if (baseline_hazard <= 0) {
    stop("baseline_hazard must be positive.", call. = FALSE)
  }
  if (censoring_max <= 0) {
    stop("censoring_max must be positive.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  numTasks <- length(thresholds)
  n_binary <- p - n_continuous

  # ---------------------------
  # Feature names
  # ---------------------------
  continuous_names <- if (n_continuous > 0) {
    paste0("cont_", seq_len(n_continuous))
  } else {
    character(0)
  }
  binary_names <- if (n_binary > 0) {
    paste0("mut_", seq_len(n_binary))
  } else {
    character(0)
  }
  feature_names <- c(continuous_names, binary_names, "treatment")
  p_total <- p + 1  # features + treatment

  # ---------------------------
  # Generate feature matrix X
  # ---------------------------
  X <- matrix(0, nrow = n, ncol = p_total)
  colnames(X) <- feature_names

  # Continuous features: standard normal
  if (n_continuous > 0) {
    X[, continuous_names] <- matrix(
      rnorm(n * n_continuous), nrow = n, ncol = n_continuous
    )
  }

  # Binary features: prevalences from Beta(2, 10)
  if (n_binary > 0) {
    prevalences <- rbeta(n_binary, shape1 = 2, shape2 = 10)
    for (j in seq_len(n_binary)) {
      X[, binary_names[j]] <- rbinom(n, size = 1, prob = prevalences[j])
    }
  }

  # Treatment: balanced 1:1
  X[, "treatment"] <- sample(rep(0:1, length.out = n))

  # ---------------------------
  # Define effect templates
  # ---------------------------
  # Each template is a vector of length numTasks, normalised to max |value| = 1
  template_library <- list(
    early      = .make_template("early", numTasks),
    late       = .make_template("late", numTasks),
    constant   = .make_template("constant", numTasks),
    switch     = .make_template("switch", numTasks),
    decreasing = .make_template("decreasing", numTasks)
  )

  template_names <- names(template_library)

  # ---------------------------
  # Assign effects to features
  # ---------------------------
  # Coefficient matrix: p_total x numTasks
  beta <- matrix(0, nrow = p_total, ncol = numTasks)
  rownames(beta) <- feature_names
  colnames(beta) <- as.character(thresholds)

  # Select which features are signals
  all_feature_idx <- seq_len(p)  # exclude treatment (handled separately)
  signal_idx <- sort(sample(all_feature_idx, n_signals))
  signal_names <- feature_names[signal_idx]
  null_idx <- setdiff(all_feature_idx, signal_idx)
  null_names <- feature_names[null_idx]

  # Assign templates cyclically with random sign
  effect_types <- character(n_signals)
  names(effect_types) <- signal_names

  for (s in seq_len(n_signals)) {
    template_choice <- template_names[((s - 1) %% length(template_names)) + 1]
    sign_choice <- sample(c(-1, 1), 1)

    beta[signal_idx[s], ] <- sign_choice * effect_strength *
      template_library[[template_choice]]

    effect_types[s] <- template_choice
  }

  # Treatment effect: constant across thresholds
  beta["treatment", ] <- treatment_effect

  # ---------------------------
  # Generate the response (mode-specific)
  # ---------------------------
  SurvTime <- NULL
  Event    <- NULL
  Y        <- NULL

  if (mode == "regression") {
    # Multi-task linear response with Gaussian noise.
    Y <- X %*% beta + matrix(rnorm(n * numTasks, sd = noise_sd),
                             nrow = n, ncol = numTasks)
    colnames(Y) <- as.character(thresholds)
  } else if (mode == "classification") {
    # Bernoulli labels in {-1, +1} from the logistic link.
    prob <- 1 / (1 + exp(-(X %*% beta)))
    Y <- matrix(ifelse(matrix(runif(n * numTasks), n, numTasks) < prob,
                       1, -1),
                nrow = n, ncol = numTasks)
    colnames(Y) <- as.character(thresholds)
  } else {
    # mode == "survival": piecewise-exponential times (historical path).
    Y <- .simulate_survival_times(
      X = X, beta = beta, n = n, numTasks = numTasks,
      thresholds = thresholds, baseline_hazard = baseline_hazard,
      censoring_max = censoring_max
    )
    SurvTime <- Y$SurvTime
    Event    <- Y$Event
    Y        <- NULL
  }

  # ---------------------------
  # Build ground truth
  # ---------------------------
  ground_truth <- list(
    coefficients    = beta,
    signal_features = signal_names,
    null_features   = null_names,
    effect_types    = effect_types,
    thresholds      = thresholds,
    baseline_hazard = baseline_hazard
  )

  # ---------------------------
  # Return S3 object
  # ---------------------------
  structure(
    list(
      mode           = mode,
      X              = X,
      Y              = Y,
      SurvTime       = SurvTime,
      Event          = Event,
      treatment      = X[, "treatment"],
      feature_names  = feature_names,
      ground_truth   = ground_truth,
      n              = n,
      p              = p,
      n_signals      = n_signals,
      n_continuous   = n_continuous,
      thresholds     = thresholds,
      seed           = seed,
      call           = cl
    ),
    class = "simulated_mtl"
  )
}


# ---------------------------------------------------------------------
# Internal: piecewise-exponential survival time generator (SIM-01)
#
# Factored out of simulate_mtl so the survival path is unchanged while
# regression/classification responses are produced inline. Returns a
# list(SurvTime, Event).
# ---------------------------------------------------------------------
.simulate_survival_times <- function(X, beta, n, numTasks, thresholds,
                                     baseline_hazard, censoring_max) {

  # Interval boundaries: [0, t1), [t1, t2), ..., [t_{K-1}, t_K), [t_K, Inf)
  interval_bounds <- c(0, thresholds)
  n_intervals <- length(interval_bounds)  # includes the tail interval

  SurvTime_true <- numeric(n)

  for (i in seq_len(n)) {
    cumulative_time <- 0
    event_occurred <- FALSE

    for (k in seq_len(n_intervals)) {
      # Coefficients for this interval
      # Use the threshold index: intervals 1..numTasks use beta[,k],
      # the tail interval (k = numTasks + 1) reuses beta[, numTasks]
      beta_k <- if (k <= numTasks) beta[, k] else beta[, numTasks]

      # Hazard for this patient in this interval
      linear_predictor <- sum(X[i, ] * beta_k)
      hazard <- baseline_hazard * exp(linear_predictor)

      # Time to event in this interval (exponential)
      time_in_interval <- rexp(1, rate = hazard)

      # Interval width (Inf for the tail)
      if (k < n_intervals) {
        interval_width <- interval_bounds[k + 1] - interval_bounds[k]
      } else {
        interval_width <- Inf
      }

      if (time_in_interval < interval_width) {
        # Event occurs in this interval
        cumulative_time <- cumulative_time + time_in_interval
        event_occurred <- TRUE
        break
      } else {
        # Survive this interval, move to next
        cumulative_time <- cumulative_time + interval_width
      }
    }

    # If no event in any interval (extremely unlikely but possible)
    if (!event_occurred) {
      cumulative_time <- cumulative_time + rexp(1, rate = baseline_hazard)
    }

    SurvTime_true[i] <- cumulative_time
  }

  # ---------------------------
  # Apply censoring
  # ---------------------------
  censor_time <- runif(n, min = 0, max = censoring_max)
  SurvTime <- pmin(SurvTime_true, censor_time)
  Event <- as.numeric(SurvTime_true <= censor_time)

  list(SurvTime = SurvTime, Event = Event)
}


# ---------------------------
# Internal: build effect templates
# ---------------------------

#' Build a normalised effect template vector
#'
#' @param type Template type: "early", "late", "constant",
#'   "increasing", "decreasing".
#' @param n_tasks Number of tasks/thresholds.
#'
#' @return Numeric vector of length \code{n_tasks}, normalised so
#'   \code{max(abs(values)) == 1}.
#'
#' @keywords internal
.make_template <- function(type, n_tasks) {

  idx <- seq_len(n_tasks)
  scaled <- (idx - 1) / max(1, n_tasks - 1)  # 0 to 1

  template <- switch(type,
                     early      = rev(scaled),
                     late       = scaled,
                     constant   = rep(1, n_tasks),
                     increasing = scaled,
                     decreasing = rev(scaled),
                     switch     = 1 - 2 * (scaled ^ 0.5),   # crosses zero early
                     stop("Unknown template type: ", type, call. = FALSE)
  )

  # Normalise to max |value| = 1
  max_abs <- max(abs(template))
  if (max_abs > 0) {
    template <- template / max_abs
  }

  return(template)
}
