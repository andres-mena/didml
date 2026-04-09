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
#' @param method Character (`"lasso"`, `"rf"`, `"nn"`, `"ols"`) or a custom
#'   function `f(Y_train, X_train, X_predict)`. Only used when
#'   `dml = TRUE`.
#' @param K Integer, number of cross-fitting folds (default 5). Only used
#'   when `dml = TRUE`.
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
#' @export
didml <- function(Y, D, G, Ti, X,
                   design = "2x2", iv = FALSE, dml = TRUE,
                   estimator = "auto", method = "lasso", K = 5L,
                   trim = "auto", trim_alpha = 0.10,
                   cluster = NULL, se_type = "analytical",
                   B = 1000L, seed = NULL, verbose = FALSE) {
  cl <- match.call()
  X <- as.matrix(X)
  N <- length(Y)
  p <- ncol(X)

  # ---- Input validation ----
  design <- match.arg(design, c("2x2", "multi"))
  if (design == "multi") {
    stop("Multi-period DID coming in a future release.", call. = FALSE)
  }

  if (iv && !dml) {
    stop("Fuzzy DID requires DML. Set dml = TRUE.", call. = FALSE)
  }

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
        K = NA_integer_,
        trim = NA_character_, trim_alpha = NA_real_,
        se_type = se_type, B = B, seed = seed
      )
    )
    class(result) <- c("didml", "didml_sharp")
    return(result)
  }

  # ---- Path 2: Sharp DID with DML (Chang 2020) ----
  if (!iv && dml) {
    estimand_nuis <- "chang"

    if (verbose) message("Step 1/4: Nuisance estimation ",
                         "(K=", K, ", method=", method, ")...")
    nuis <- didml_nuisance(Y, D, G, Ti, X,
                            method = method, K = K,
                            iv = FALSE, seed = seed)

    if (verbose) message("Step 2/4: Propensity trimming (", trim, ")...")
    # For sharp DID, no DID_D needed for trimming
    trim_info <- didml_trim(nuis$pG_raw,
                             DID_D = NULL,
                             method = trim,
                             alpha_fixed = trim_alpha)

    # Apply trimming to propensity
    nuis$pG_raw_original <- nuis$pG_raw
    pG_trimmed <- trim_info$pG_trimmed
    pT <- nuis$pT
    nuis$pG_raw <- pG_trimmed
    nuis$pi_11 <- pG_trimmed * pT
    nuis$pi_10 <- pG_trimmed * (1 - pT)
    nuis$pi_01 <- (1 - pG_trimmed) * pT
    nuis$pi_00 <- (1 - pG_trimmed) * (1 - pT)

    if (verbose) message("Step 3/4: Computing Chang (2020) estimator...")
    W <- data.frame(Y = Y, D = D, G = G, Ti = Ti)
    chang_res <- didml_chang(W, nuis)

    chang_inf <- didml_inference(chang_res,
                                  cluster = cluster,
                                  se_type = se_type,
                                  B = B, seed = seed)

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

    if (verbose) message("Step 4/4: Done.")

    result <- list(
      estimates = estimates,
      wald = NULL, tc = NULL,
      chang = chang_res,
      drdid = NULL,
      nuisance = nuis,
      trim_info = trim_info,
      call = cl,
      N = N, p = p,
      design = design, iv = iv, dml = dml,
      settings = list(
        estimator = "chang", method = method, K = K,
        trim = trim, trim_alpha = trim_alpha,
        se_type = se_type, B = B, seed = seed
      )
    )
    class(result) <- c("didml", "didml_sharp")
    return(result)
  }

  # ---- Path 3: Fuzzy DID with DML (Wald / TC) ----
  estimand <- match.arg(estimator, c("wald", "tc", "both"))

  if (verbose) message("Step 1/4: Nuisance estimation ",
                       "(K=", K, ", method=", method, ")...")
  nuis <- didml_nuisance(Y, D, G, Ti, X,
                          method = method, K = K,
                          iv = TRUE, estimand = estimand,
                          seed = seed)

  if (verbose) message("Step 2/4: Propensity trimming (", trim, ")...")
  W <- data.frame(Y = Y, D = D, G = G, Ti = Ti)

  DID_D_hat <- D - nuis$m_D_10 - nuis$m_D_01 + nuis$m_D_00
  trim_info <- didml_trim(nuis$pG_raw, DID_D = DID_D_hat,
                           method = trim, alpha_fixed = trim_alpha)

  # Apply trimming to propensity — preserve raw for diagnostics
  pT <- mean(Ti)
  nuis$pG_raw_original <- nuis$pG_raw
  pG_trimmed <- trim_info$pG_trimmed
  nuis$pi_11 <- pG_trimmed * pT
  nuis$pi_10 <- pG_trimmed * (1 - pT)
  nuis$pi_01 <- (1 - pG_trimmed) * pT
  nuis$pi_00 <- (1 - pG_trimmed) * (1 - pT)

  if (verbose) message("Step 3/4: Computing estimators...")
  wald_res <- tc_res <- NULL
  estimates <- data.frame()

  if (estimand %in% c("wald", "both")) {
    wald_res <- didml_wald(W, nuis)
    wald_inf <- didml_inference(wald_res, cluster = cluster,
                                se_type = se_type, B = B, seed = seed)
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

  if (estimand %in% c("tc", "both")) {
    tc_res <- didml_tc(W, nuis)
    tc_inf <- didml_inference(tc_res, cluster = cluster,
                              se_type = se_type, B = B, seed = seed)
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

  if (verbose) message("Step 4/4: Done.")

  result <- list(
    estimates = estimates,
    wald = wald_res, tc = tc_res,
    chang = NULL,
    drdid = NULL,
    nuisance = nuis,
    trim_info = trim_info,
    call = cl,
    N = N, p = p,
    design = design, iv = iv, dml = dml,
    settings = list(
      estimator = estimand, method = method, K = K,
      trim = trim, trim_alpha = trim_alpha,
      se_type = se_type, B = B, seed = seed
    )
  )
  class(result) <- c("didml", "didml_fuzzy")
  result
}
