#' Optimal Propensity Score Trimming
#'
#' Implements the data-driven one-sided trimming rule from Corollary 2.1
#' and Remark 2 of Mena (2026). Minimizes the estimated variance objective
#' R(alpha) over a grid of trimming thresholds.
#'
#' @param pG_raw Numeric vector of cross-fitted propensity scores Pr(G=1|X).
#' @param DID_D Numeric vector of first-stage DID predictions per observation.
#' @param method Character: `"auto"` (grid search, default), `"fixed"`, or
#'   `"none"`.
#' @param alpha_fixed Numeric trimming threshold for `method = "fixed"`
#'   (default 0.10, symmetric).
#' @param alpha_grid Numeric vector of candidate alpha values for grid search.
#'
#' @return A list with:
#'   \describe{
#'     \item{pG_trimmed}{Trimmed propensity scores.}
#'     \item{alpha_hat}{Selected trimming threshold (0 if no trimming).}
#'     \item{keep}{Logical vector indicating which observations survive trimming.}
#'   }
#'
#' @export
didml_trim <- function(pG_raw,
                        DID_D = NULL,
                        method = "auto",
                        alpha_fixed = 0.10,
                        alpha_grid = seq(0, 0.50, by = 0.01)) {
  N <- length(pG_raw)

  if (method == "none") {
    return(list(pG_trimmed = pG_raw, alpha_hat = 0, keep = rep(TRUE, N)))
  }

  if (method == "fixed") {
    pG_trimmed <- pmin(pG_raw, 1 - alpha_fixed)
    pG_trimmed <- pmax(pG_trimmed, alpha_fixed)
    keep <- pG_raw >= alpha_fixed & pG_raw <= (1 - alpha_fixed)
    return(list(pG_trimmed = pG_trimmed, alpha_hat = alpha_fixed, keep = keep))
  }

  # method == "auto": data-driven one-sided trimming (Corollary 2.1, Remark 2)
  # Estimate g(e) = E[DID_D | e(X) = e] via local linear regression
  if (is.null(DID_D)) {
    # Without first-stage predictions, use Corollary 2.2 (constant g)
    g_hat <- rep(1, N)
  } else {
    g_fit <- tryCatch(
      stats::loess(DID_D ~ pG_raw, span = 0.3, degree = 1),
      error = function(e) NULL
    )
    if (!is.null(g_fit)) {
      g_hat <- pmax(stats::predict(g_fit), 0.001)
      g_hat[is.na(g_hat)] <- mean(DID_D, na.rm = TRUE)
    } else {
      g_hat <- rep(mean(DID_D, na.rm = TRUE), N)
    }
  }

  # Grid search: minimize R_hat(alpha)
  R_vals <- numeric(length(alpha_grid))
  for (i in seq_along(alpha_grid)) {
    a <- alpha_grid[i]
    sel <- pG_raw <= (1 - a)
    if (sum(sel) < 10) { R_vals[i] <- Inf; next }
    num <- sum(pG_raw[sel] / (1 - pG_raw[sel]))
    den <- sum(pG_raw[sel] * g_hat[sel])^2
    R_vals[i] <- if (den > 1e-10) num / den else Inf
  }

  alpha_hat <- alpha_grid[which.min(R_vals)]
  pG_trimmed <- pmin(pG_raw, 1 - alpha_hat)
  keep <- pG_raw <= (1 - alpha_hat)

  list(pG_trimmed = pG_trimmed, alpha_hat = alpha_hat, keep = keep)
}
