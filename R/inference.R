#' Inference for didml Results
#'
#' Computes standard errors, confidence intervals, and p-values for
#' DML-Wald, DML-TC, and DML-Chang estimates. Supports analytical (i.i.d.),
#' clustered, and multiplier bootstrap variance estimation.
#'
#' @param result A list from \code{didml_wald()}, \code{didml_tc()}, or
#'   \code{didml_chang()}, must contain `estimate`, `scores`, `den_i`,
#'   and `variance`.
#' @param cluster Optional numeric or factor vector of cluster identifiers
#'   (same length as the data). When provided, computes cluster-robust SEs.
#' @param se_type Character: `"analytical"` (default), `"cluster"` (requires
#'   `cluster`), or `"bootstrap"`.
#' @param B Integer, number of bootstrap replications (default 1000).
#' @param alpha Significance level for confidence intervals (default 0.05).
#' @param seed Random seed for bootstrap (default `NULL`).
#'
#' @return A list with `estimate`, `se`, `ci_lower`, `ci_upper`, `p_value`,
#'   `variance`, `se_type`, `alpha`.
#'
#' @export
didml_inference <- function(result, cluster = NULL, se_type = "analytical",
                             B = 1000L, alpha = 0.05, seed = NULL) {
  se_type <- match.arg(se_type, c("analytical", "cluster", "bootstrap"))

  # Validate result object
  if (is.null(result$den_i)) stop("`result` must contain `den_i`.", call. = FALSE)
  if (is.null(result$scores)) stop("`result` must contain `scores`.", call. = FALSE)

  theta_hat <- result$estimate
  if (!is.finite(theta_hat)) {
    return(list(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_,
                ci_upper = NA_real_, p_value = NA_real_, variance = NA_real_,
                se_type = se_type, alpha = alpha))
  }

  psi_i <- result$scores
  N <- length(psi_i)
  G_hat <- -mean(result$den_i, na.rm = TRUE)
  z_alpha <- stats::qnorm(1 - alpha / 2)

  # Resolve cluster vs se_type conflict

  if (!is.null(cluster) && se_type != "cluster") {
    warning("`cluster` provided; switching to cluster-robust SEs.", call. = FALSE)
    se_type <- "cluster"
  }

  if (se_type == "cluster") {
    if (is.null(cluster))
      stop("`cluster` must be provided for se_type = 'cluster'.", call. = FALSE)
    if (length(cluster) != N)
      stop("`cluster` must have the same length as the data.", call. = FALSE)

    # Cluster-robust variance with S/(S-1) finite-sample correction
    cluster_sums <- tapply(psi_i, cluster, function(x) sum(x, na.rm = TRUE))
    S <- length(cluster_sums)
    if (S < 10) {
      warning("Only ", S, " clusters. Cluster-robust SEs may be unreliable with fewer than 10.",
              call. = FALSE)
    }
    V_hat <- (S / (S - 1)) * sum(cluster_sums^2) / (N^2 * G_hat^2)
    se_hat <- sqrt(V_hat)

  } else if (se_type == "bootstrap") {
    # Score-based multiplier bootstrap (Chernozhukov et al. 2018, Section 4.2):
    # theta_b = theta_hat + (1/N) * sum(w_b * psi_i / G_hat)
    # where w_b are i.i.d. Rademacher weights.
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
        get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        if (!is.null(old_seed))
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        else
          rm(".Random.seed", envir = .GlobalEnv)
      })
      set.seed(seed)
    }
    ok <- is.finite(psi_i)
    psi_ok <- psi_i[ok]
    N_ok <- length(psi_ok)
    boot_theta <- numeric(B)
    for (b in seq_len(B)) {
      w_b <- sample(c(-1, 1), N_ok, replace = TRUE)
      boot_theta[b] <- theta_hat + sum(w_b * psi_ok) / (N_ok * G_hat)
    }
    se_hat <- stats::sd(boot_theta)
    V_hat <- se_hat^2

  } else {
    # Analytical (i.i.d.)
    V_hat <- result$variance
    N_ok <- sum(is.finite(psi_i))
    se_hat <- sqrt(V_hat / N_ok)
  }

  ci_lower <- theta_hat - z_alpha * se_hat
  ci_upper <- theta_hat + z_alpha * se_hat
  p_value <- 2 * stats::pnorm(-abs(theta_hat / se_hat))

  list(
    estimate = theta_hat,
    se = se_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    variance = V_hat,
    se_type = se_type,
    alpha = alpha
  )
}
