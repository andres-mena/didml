#' Internal DRDID Wrapper for Sharp DID without DML
#'
#' Wraps \code{DRDID::drdid_rc()} for repeated cross-section doubly robust
#' DID estimation (Sant'Anna & Zhao 2020). Critical mapping: DRDID's
#' \code{D} argument = our \code{G} (group), DRDID's \code{post} = our
#' \code{Ti} (time period).
#'
#' @param Y Numeric vector of outcomes.
#' @param D Numeric vector of treatment take-up.
#' @param G Binary vector of group assignment (0/1).
#' @param Ti Binary vector of time period (0/1).
#' @param X Numeric matrix or data.frame of covariates.
#' @param cluster Optional cluster identifiers (not used by DRDID directly,
#'   reserved for future use).
#' @param se_type Character: \code{"analytical"} (default) or
#'   \code{"bootstrap"}.
#' @param B Integer, bootstrap replications (default 1000).
#' @param seed Random seed (default \code{NULL}).
#' @param verbose Logical, print progress (default \code{FALSE}).
#'
#' @return A list with \code{estimate}, \code{se}, \code{ci_lower},
#'   \code{ci_upper}, \code{scores}, \code{p_value}, \code{drdid_raw}.
#'
#' @noRd
.didml_drdid <- function(Y, D, G, Ti, X, cluster = NULL,
                          se_type = "analytical", B = 1000L,
                          seed = NULL, verbose = FALSE) {
  if (!requireNamespace("DRDID", quietly = TRUE))
    stop("Install 'DRDID': install.packages('DRDID')", call. = FALSE)

  # Warn if D != G*Ti for many obs (fuzzy compliance)
  compliance_rate <- mean(D == (G * Ti), na.rm = TRUE)
  if (compliance_rate < 0.95)
    warning("D differs from G*Ti for ",
            round(100 * (1 - compliance_rate), 1),
            "% of obs. DRDID estimates the sharp ATT on group G, ",
            "ignoring treatment take-up D. ",
            "For fuzzy designs, use iv = TRUE.",
            call. = FALSE)

  # DRDID::drdid_rc expects covariates WITH intercept column
  # (it does not add one internally unlike the high-level drdid())
  X_mat <- cbind(1, as.matrix(X))
  use_boot <- (se_type == "bootstrap")

  drdid_out <- DRDID::drdid_rc(
    y = Y, post = Ti, D = G,
    covariates = X_mat,
    boot = use_boot,
    nboot = if (use_boot) B else NULL,
    inffunc = TRUE
  )

  est <- drdid_out$ATT
  se_val <- drdid_out$se

  list(
    estimate = est,
    se = se_val,
    ci_lower = drdid_out$lci,
    ci_upper = drdid_out$uci,
    scores = if (!is.null(drdid_out$att.inf.func))
      drdid_out$att.inf.func else NULL,
    p_value = 2 * stats::pnorm(-abs(est / se_val)),
    drdid_raw = drdid_out
  )
}
