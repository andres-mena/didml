#' Difference-in-Differences with Machine Learning
#'
#' Unified estimation for difference-in-differences designs with machine
#' learning nuisance estimation. Supports sharp DID (Sant'Anna & Zhao 2020
#' via DRDID, and DML-DID of Chang 2020), fuzzy DID (DML-Wald and DML-TC
#' of Mena 2026), and extensible multi-period designs.
#'
#' @param Y Numeric vector of outcomes.
#' @param D Numeric vector of treatment (binary or continuous).
#' @param G Binary vector of group assignment (0/1).
#' @param Ti Binary vector of time period (0/1).
#' @param X Numeric matrix or data.frame of covariates.
#' @param design Character: `"2x2"` (default) for two-group, two-period
#'   designs, or `"multi"` for multi-period (not yet implemented).
#' @param iv Logical. If `TRUE`, use fuzzy DID estimators (DML-Wald/TC)
#'   that account for imperfect compliance. If `FALSE` (default), use
#'   sharp DID estimators.
#' @param dml Logical. If `TRUE` (default), use double/debiased machine
#'   learning for nuisance estimation. If `FALSE`, use parametric DRDID
#'   (requires \code{iv = FALSE}).
#' @param estimator Character: `"auto"` (default), `"wald"`, `"tc"`,
#'   `"both"`, `"chang"`, or `"drdid"`. When `"auto"`:
#'   \itemize{
#'     \item \code{iv = TRUE}: resolves to `"both"` (Wald + TC).
#'     \item \code{iv = FALSE, dml = TRUE}: resolves to `"chang"`.
#'     \item \code{iv = FALSE, dml = FALSE}: resolves to `"drdid"`.
#'   }
#' @param method Character string, character vector, or custom function.
#'   A single string (`"lasso"`, `"rf"`, `"nn"`, `"ols"`) or a function
#'   `f(Y_train, X_train, X_predict)` uses one learner. A character vector
#'   (e.g. `c("ols", "lasso", "rf")`) activates stacking, combining
#'   predictions via \code{ensemble_type}. Only used when `dml = TRUE`.
#' @param method_outcome Character string or vector overriding \code{method}
#'   for outcome/treatment regressions (m_Y, m_D, ell_20). Default
#'   \code{NULL} uses \code{method}.
#' @param method_propensity Character string or vector overriding \code{method}
#'   for propensity score estimation. Default \code{NULL} uses \code{method}.
#' @param ensemble_type Character: `"average"` (default when multiple methods),
#'   `"nnls1"` (non-negative least squares, weights sum to 1), or
#'   `"singlebest"` (pick learner with lowest CV MSE). Ignored for single
#'   learner.
#' @param K Integer, number of cross-fitting folds (default 5). Only used
#'   when `dml = TRUE`.
#' @param reps Integer, number of cross-fitting repetitions (default 1).
#'   When \code{reps > 1}, the full pipeline (nuisance, trim, score,
#'   inference) is run \code{reps} times with independent fold assignments.
#'   Estimates are aggregated using the median; SEs incorporate both
#'   within-rep variance and across-rep variability (DoubleML convention).
#' @param trim Character: `"auto"` (data-driven, default), `"fixed"`,
#'   or `"none"`. Only used when `dml = TRUE`.
#' @param trim_alpha Numeric, fixed trimming level when `trim = "fixed"`.
#' @param cluster Optional vector of cluster identifiers for cluster-robust
#'   SEs.
#' @param se_type Character: `"analytical"` (default), `"cluster"`, or
#'   `"bootstrap"`.
#' @param B Integer, bootstrap replications (default 1000).
#' @param seed Random seed for cross-fitting folds and bootstrap
#'   (default `NULL`).
#' @param verbose Logical, print progress messages (default `FALSE`).
#'
#' @return An object of class `"didml"` (with subclass `"didml_fuzzy"` or
#'   `"didml_sharp"`) containing:
#'   \describe{
#'     \item{estimates}{Data frame with columns: estimand, estimate, se,
#'       ci_lower, ci_upper, p_value.}
#'     \item{wald}{Full DML-Wald results (if requested, else NULL).}
#'     \item{tc}{Full DML-TC results (if requested, else NULL).}
#'     \item{chang}{Full DML-Chang results (if requested, else NULL).}
#'     \item{drdid}{Full DRDID results (if requested, else NULL).}
#'     \item{nuisance}{Cross-fitted nuisance predictions (if DML).}
#'     \item{trim_info}{Trimming diagnostics (if DML).}
#'     \item{call}{The matched call.}
#'     \item{N}{Sample size.}
#'     \item{p}{Number of covariates.}
#'     \item{design}{The design type.}
#'     \item{iv}{Whether fuzzy DID was used.}
#'     \item{dml}{Whether DML was used.}
#'     \item{settings}{List of estimation settings.}
#'     \item{ensemble_weights}{Ensemble weights when stacking (NULL otherwise).}
#'     \item{reps_detail}{Per-rep results when \code{reps > 1} (NULL otherwise).}
#'   }
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' N <- 500
#' X <- matrix(rnorm(N * 5), ncol = 5)
#' G <- rbinom(N, 1, plogis(0.5 * X[,1]))
#' Ti <- rbinom(N, 1, 0.5)
#' D <- rbinom(N, 1, plogis(X[,1] + G + Ti * G))
#' Y <- 1 + X[,1] + 0.3 * G * Ti + D * (G * Ti) + rnorm(N)
#'
#' # Fuzzy DID (DML-Wald)
#' fit_fuzzy <- didml(Y, D, G, Ti, X, iv = TRUE,
#'                    estimator = "wald", method = "ols", K = 2)
#' print(fit_fuzzy)
#'
#' # Sharp DID (DML-Chang)
#' fit_sharp <- didml(Y, D, G, Ti, X, iv = FALSE, dml = TRUE,
#'                    method = "ols", K = 2)
#' print(fit_sharp)
#'
#' # Stacking with multiple learners
#' fit_stack <- didml(Y, D, G, Ti, X, iv = FALSE, dml = TRUE,
#'                    method = c("ols", "lasso"), K = 2)
#'
#' # Cross-fitting repetitions
#' fit_reps <- didml(Y, D, G, Ti, X, iv = FALSE, dml = TRUE,
#'                   method = "ols", K = 2, reps = 3)
#'
#' # Separate learners for outcome and propensity
#' fit_sep <- didml(Y, D, G, Ti, X, iv = FALSE, dml = TRUE,
#'                  method = "ols", method_propensity = "lasso", K = 2)
#' }
#'
#' @references
#' Mena, A. (2026). Double Debiased Machine Learning for
#' Difference-in-Differences under Imperfect Compliance. Working Paper,
#' Brown University.
#'
#' De Chaisemartin, C. and D'Haultfoeuille, X. (2018). Fuzzy
#' Differences-in-Differences. \emph{Review of Economic Studies}, 85(2),
#' 999-1028.
#'
#' Chang, N.-C. (2020). Double/debiased machine learning for
#' difference-in-differences treatment effects. \emph{The Econometrics
#' Journal}, 23(2), 177-191.
#'
#' Sant'Anna, P. H. C. and Zhao, J. (2020). Doubly robust
#' difference-in-differences estimators. \emph{Journal of Econometrics},
#' 219(1), 101-122.
#'
#' Ahrens, A., Hansen, C. B., Schaffer, M. E. and Wiemann, T. (2024).
#' ddml: Double/debiased machine learning in Stata. \emph{The Stata
#' Journal}, 24(1), 3-45.
#'
#' @export
didml <- function(Y, D, G, Ti, X,
                   design = "2x2", iv = FALSE, dml = TRUE,
                   estimator = "auto",
                   method = "lasso",
                   method_outcome = NULL,
                   method_propensity = NULL,
                   ensemble_type = "average",
                   K = 5L, reps = 1L,
                   trim = "auto", trim_alpha = 0.10,
                   cluster = NULL, se_type = "analytical",
                   B = 1000L, seed = NULL, verbose = FALSE) {
  cl <- match.call()
  X <- as.matrix(X)
  N <- length(Y)
  p <- ncol(X)
  reps <- as.integer(reps)

  # ---- Input validation ----
  design <- match.arg(design, c("2x2", "multi"))
  if (design == "multi") {
    stop("Multi-period DID coming in a future release.", call. = FALSE)
  }

  if (iv && !dml) {
    stop("Fuzzy DID requires DML. Set dml = TRUE.", call. = FALSE)
  }

  if (reps < 1L) stop("`reps` must be a positive integer.", call. = FALSE)

  .validate_inputs(Y, D, G, Ti, X)

  # ---- Resolve estimator ----
  if (estimator == "auto") {
    if (iv) {
      estimator <- "both"
    } else if (dml) {
      estimator <- "chang"
    } else {
      estimator <- "drdid"
    }
  }

  # ---- Path 1: DRDID (sharp, no DML) ----
  if (!iv && !dml) {
    if (verbose) message("Estimating sharp DID via DRDID...")
    drdid_res <- .didml_drdid(Y, D, G, Ti, X,
                               cluster = cluster,
                               se_type = se_type,
                               B = B, seed = seed,
                               verbose = verbose)

    estimates <- data.frame(
      estimand = "DRDID",
      estimate = drdid_res$estimate,
      se = drdid_res$se,
      ci_lower = drdid_res$ci_lower,
      ci_upper = drdid_res$ci_upper,
      p_value = drdid_res$p_value,
      stringsAsFactors = FALSE
    )

    result <- list(
      estimates = estimates,
      wald = NULL, tc = NULL, chang = NULL,
      drdid = drdid_res,
      nuisance = NULL,
      trim_info = NULL,
      call = cl,
      N = N, p = p,
      design = design, iv = iv, dml = dml,
      settings = list(
        estimator = "drdid", method = NA_character_,
        K = NA_integer_, reps = 1L,
        trim = NA_character_, trim_alpha = NA_real_,
        se_type = se_type, B = B, seed = seed
      )
    )
    class(result) <- c("didml", "didml_sharp")
    return(result)
  }

  # ---- Helper: run a single DML rep ----
  .run_single_rep <- function(rep_seed) {
    nuis_args <- list(
      Y = Y, D = D, G = G, Ti = Ti, X = X,
      method = method,
      method_outcome = method_outcome,
      method_propensity = method_propensity,
      ensemble_type = ensemble_type,
      K = K, seed = rep_seed
    )

    # ---- Sharp DID with DML (Chang 2020) ----
    if (!iv && dml) {
      nuis_args$iv <- FALSE
      nuis <- do.call(didml_nuisance, nuis_args)

      trim_info <- didml_trim(nuis$pG_raw, DID_D = NULL,
                               method = trim, alpha_fixed = trim_alpha)

      nuis$pG_raw_original <- nuis$pG_raw
      pG_trimmed <- trim_info$pG_trimmed
      pT <- nuis$pT
      nuis$pG_raw <- pG_trimmed
      nuis$pi_11 <- pG_trimmed * pT
      nuis$pi_10 <- pG_trimmed * (1 - pT)
      nuis$pi_01 <- (1 - pG_trimmed) * pT
      nuis$pi_00 <- (1 - pG_trimmed) * (1 - pT)

      W <- data.frame(Y = Y, D = D, G = G, Ti = Ti)
      chang_res <- didml_chang(W, nuis)
      chang_inf <- didml_inference(chang_res, cluster = cluster,
                                    se_type = se_type, B = B, seed = rep_seed)

      estimates <- data.frame(
        estimand = "DML-Chang",
        estimate = chang_inf$estimate,
        se = chang_inf$se,
        ci_lower = chang_inf$ci_lower,
        ci_upper = chang_inf$ci_upper,
        p_value = chang_inf$p_value,
        stringsAsFactors = FALSE
      )

      chang_res <- c(chang_res,
                     chang_inf[c("se", "ci_lower", "ci_upper",
                                 "p_value", "se_type")])

      return(list(
        estimates = estimates,
        wald = NULL, tc = NULL,
        chang = chang_res, drdid = NULL,
        nuisance = nuis, trim_info = trim_info,
        subclass = "didml_sharp",
        estimand_resolved = "chang"
      ))
    }

    # ---- Fuzzy DID with DML (Wald / TC) ----
    estimand_f <- match.arg(estimator, c("wald", "tc", "both"))
    nuis_args$iv <- TRUE
    nuis_args$estimand <- estimand_f
    nuis <- do.call(didml_nuisance, nuis_args)

    W <- data.frame(Y = Y, D = D, G = G, Ti = Ti)
    DID_D_hat <- D - nuis$m_D_10 - nuis$m_D_01 + nuis$m_D_00
    trim_info <- didml_trim(nuis$pG_raw, DID_D = DID_D_hat,
                             method = trim, alpha_fixed = trim_alpha)

    pT_val <- mean(Ti)
    nuis$pG_raw_original <- nuis$pG_raw
    pG_trimmed <- trim_info$pG_trimmed
    nuis$pi_11 <- pG_trimmed * pT_val
    nuis$pi_10 <- pG_trimmed * (1 - pT_val)
    nuis$pi_01 <- (1 - pG_trimmed) * pT_val
    nuis$pi_00 <- (1 - pG_trimmed) * (1 - pT_val)

    wald_res <- tc_res <- NULL
    estimates <- data.frame()

    if (estimand_f %in% c("wald", "both")) {
      wald_res <- didml_wald(W, nuis)
      wald_inf <- didml_inference(wald_res, cluster = cluster,
                                  se_type = se_type, B = B, seed = rep_seed)
      estimates <- rbind(estimates, data.frame(
        estimand = "DML-Wald",
        estimate = wald_inf$estimate,
        se = wald_inf$se,
        ci_lower = wald_inf$ci_lower,
        ci_upper = wald_inf$ci_upper,
        p_value = wald_inf$p_value,
        stringsAsFactors = FALSE
      ))
      wald_res <- c(wald_res,
                    wald_inf[c("se", "ci_lower", "ci_upper",
                               "p_value", "se_type")])
    }

    if (estimand_f %in% c("tc", "both")) {
      tc_res <- didml_tc(W, nuis)
      tc_inf <- didml_inference(tc_res, cluster = cluster,
                                se_type = se_type, B = B, seed = rep_seed)
      estimates <- rbind(estimates, data.frame(
        estimand = "DML-TC",
        estimate = tc_inf$estimate,
        se = tc_inf$se,
        ci_lower = tc_inf$ci_lower,
        ci_upper = tc_inf$ci_upper,
        p_value = tc_inf$p_value,
        stringsAsFactors = FALSE
      ))
      tc_res <- c(tc_res,
                  tc_inf[c("se", "ci_lower", "ci_upper",
                           "p_value", "se_type")])
    }

    list(
      estimates = estimates,
      wald = wald_res, tc = tc_res,
      chang = NULL, drdid = NULL,
      nuisance = nuis, trim_info = trim_info,
      subclass = "didml_fuzzy",
      estimand_resolved = estimand_f
    )
  }

  # ---- Run reps ----
  if (reps == 1L) {
    # Single rep — backward compatible path
    rep_seed <- seed
    if (verbose) message("Step 1/4: Nuisance estimation ",
                         "(K=", K, ", method=",
                         paste(method, collapse = "+"), ")...")
    single <- .run_single_rep(rep_seed)

    if (verbose) message("Step 4/4: Done.")

    # Format method for settings
    method_display <- if (is.character(method) && length(method) > 1) {
      paste(method, collapse = "+")
    } else if (is.character(method)) {
      method
    } else {
      "custom"
    }

    result <- list(
      estimates = single$estimates,
      wald = single$wald, tc = single$tc,
      chang = single$chang, drdid = single$drdid,
      nuisance = single$nuisance,
      trim_info = single$trim_info,
      call = cl,
      N = N, p = p,
      design = design, iv = iv, dml = dml,
      settings = list(
        estimator = single$estimand_resolved,
        method = method_display,
        method_outcome = method_outcome,
        method_propensity = method_propensity,
        ensemble_type = if (is.character(method) && length(method) > 1) ensemble_type else NULL,
        K = K, reps = 1L,
        trim = trim, trim_alpha = trim_alpha,
        se_type = se_type, B = B, seed = seed
      ),
      ensemble_weights = single$nuisance$ensemble_weights,
      reps_detail = NULL
    )
    class(result) <- c("didml", single$subclass)
    return(result)

  } else {
    # ---- Multiple reps ----
    rep_results <- vector("list", reps)
    for (r in seq_len(reps)) {
      rep_seed <- if (!is.null(seed)) seed + r else NULL
      if (verbose) message("Rep ", r, "/", reps, ": running pipeline...")
      rep_results[[r]] <- .run_single_rep(rep_seed)
    }

    # Aggregate across reps
    # Collect per-estimand results
    all_estimands <- unique(rep_results[[1]]$estimates$estimand)
    agg_rows <- list()

    for (est_name in all_estimands) {
      coefs <- vapply(rep_results, function(rr) {
        row <- rr$estimates[rr$estimates$estimand == est_name, ]
        if (nrow(row) == 0) return(NA_real_)
        row$estimate
      }, numeric(1))

      ses <- vapply(rep_results, function(rr) {
        row <- rr$estimates[rr$estimates$estimand == est_name, ]
        if (nrow(row) == 0) return(NA_real_)
        row$se
      }, numeric(1))

      # Aggregate: median of coefficients
      med_coef <- stats::median(coefs, na.rm = TRUE)

      # Aggregate SE: median of sqrt(se^2 + (coef - median_coef)^2) (DoubleML convention)
      agg_se <- stats::median(sqrt(ses^2 + (coefs - med_coef)^2), na.rm = TRUE)

      z_alpha <- stats::qnorm(0.975)
      agg_rows[[est_name]] <- data.frame(
        estimand = est_name,
        estimate = med_coef,
        se = agg_se,
        ci_lower = med_coef - z_alpha * agg_se,
        ci_upper = med_coef + z_alpha * agg_se,
        p_value = 2 * stats::pnorm(-abs(med_coef / agg_se)),
        stringsAsFactors = FALSE
      )
    }

    estimates <- do.call(rbind, agg_rows)
    rownames(estimates) <- NULL

    # Use last rep's nuisance/trim for diagnostics
    last <- rep_results[[reps]]

    method_display <- if (is.character(method) && length(method) > 1) {
      paste(method, collapse = "+")
    } else if (is.character(method)) {
      method
    } else {
      "custom"
    }

    # Build per-rep detail
    reps_detail <- lapply(seq_len(reps), function(r) {
      rep_results[[r]]$estimates
    })

    result <- list(
      estimates = estimates,
      wald = last$wald, tc = last$tc,
      chang = last$chang, drdid = last$drdid,
      nuisance = last$nuisance,
      trim_info = last$trim_info,
      call = cl,
      N = N, p = p,
      design = design, iv = iv, dml = dml,
      settings = list(
        estimator = last$estimand_resolved,
        method = method_display,
        method_outcome = method_outcome,
        method_propensity = method_propensity,
        ensemble_type = if (is.character(method) && length(method) > 1) ensemble_type else NULL,
        K = K, reps = reps,
        trim = trim, trim_alpha = trim_alpha,
        se_type = se_type, B = B, seed = seed
      ),
      ensemble_weights = last$nuisance$ensemble_weights,
      reps_detail = reps_detail
    )
    class(result) <- c("didml", last$subclass)

    if (verbose) message("Aggregated ", reps, " reps. Done.")
    return(result)
  }
}
