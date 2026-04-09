#' Compute DML-Wald Scores
#'
#' Computes the locally robust augmented score for the Wald-DID estimand.
#' Sign convention: (-, -, +) for correction terms w_10, w_01, w_00,
#' matching DID = Y - m_10 - m_01 + m_00.
#' Verified against CD'H (2018) supplement p.45 (r = 1.0000).
#'
#' @param W Data frame with columns Y, D, G, Ti.
#' @param pred List of cross-fitted nuisance predictions from \code{didml_nuisance()}.
#'
#' @return List with `estimate`, `variance`, `scores`, `num_i`, `den_i`,
#'   `wald_naive`. Returns NA estimate with warning if first stage is too weak.
#'
#' @export
didml_wald <- function(W, pred) {
  N <- nrow(W)
  Y <- W$Y; D <- W$D; G <- W$G; Ti <- W$Ti
  GT <- G * Ti
  p_11 <- mean(GT)

  if (p_11 < 1e-10) {
    warning("No observations in treated group at post-period (G=1, T=1).", call. = FALSE)
    return(list(estimate = NA_real_, variance = NA_real_, scores = rep(NA_real_, N),
                num_i = rep(NA_real_, N), den_i = rep(NA_real_, N), wald_naive = NA_real_))
  }

  # DID of conditional expectations
  DID_Y <- Y - pred$m_Y_10 - pred$m_Y_01 + pred$m_Y_00
  DID_D <- D - pred$m_D_10 - pred$m_D_01 + pred$m_D_00

  wald_naive <- sum(GT * DID_Y, na.rm = TRUE) / sum(GT * DID_D, na.rm = TRUE)

  # Riesz representer weights (CD'H 2018, Supplement p.45)
  w_10 <- G * (1 - Ti) * pred$pi_11 / (pred$pi_10 * p_11)
  w_01 <- (1 - G) * Ti * pred$pi_11 / (pred$pi_01 * p_11)
  w_00 <- (1 - G) * (1 - Ti) * pred$pi_11 / (pred$pi_00 * p_11)

  # Guard non-finite weights — set to NA, not zero
  for (nm in c("w_10", "w_01", "w_00")) {
    w <- get(nm)
    if (any(!is.finite(w))) {
      warning("Non-finite weights in ", nm, " (near-zero propensity). ",
              "Affected obs set to NA.", call. = FALSE)
      w[!is.finite(w)] <- NA_real_
      assign(nm, w)
    }
  }

  # Residuals
  r_Y_10 <- Y - pred$m_Y_10; r_Y_01 <- Y - pred$m_Y_01; r_Y_00 <- Y - pred$m_Y_00
  r_D_10 <- D - pred$m_D_10; r_D_01 <- D - pred$m_D_01; r_D_00 <- D - pred$m_D_00

  # Closed-form GMM: numerator and denominator
  num_i <- (GT / p_11) * DID_Y - w_10 * r_Y_10 - w_01 * r_Y_01 + w_00 * r_Y_00
  den_i <- (GT / p_11) * DID_D - w_10 * r_D_10 - w_01 * r_D_01 + w_00 * r_D_00

  sum_den <- sum(den_i, na.rm = TRUE)
  if (abs(sum_den) < 1e-10) {
    warning("Weak first stage: denominator near zero. Wald ratio unstable. Returning NA.", call. = FALSE)
    return(list(estimate = NA_real_, variance = NA_real_, scores = rep(NA_real_, N),
                num_i = num_i, den_i = den_i, wald_naive = wald_naive))
  }

  theta_hat <- sum(num_i, na.rm = TRUE) / sum_den

  # Scores at theta_hat
  psi_i <- (GT / p_11) * (DID_Y - DID_D * theta_hat) -
    w_10 * (r_Y_10 - r_D_10 * theta_hat) -
    w_01 * (r_Y_01 - r_D_01 * theta_hat) +
    w_00 * (r_Y_00 - r_D_00 * theta_hat)

  # Variance (Theorem 1)
  G_hat <- -mean(den_i, na.rm = TRUE)
  ok <- is.finite(psi_i)
  V_hat <- mean(psi_i[ok]^2) / (G_hat^2)

  list(
    estimate = theta_hat,
    variance = V_hat,
    scores = psi_i,
    num_i = num_i,
    den_i = den_i,
    wald_naive = wald_naive
  )
}


#' Compute DML-TC Scores
#'
#' Computes the locally robust augmented score for the Time-Corrected estimand.
#' Equivalent to CD'H (2018) supplement p.45 up to normalization.
#'
#' @param W Data frame with columns Y, D, G, Ti.
#' @param pred List of cross-fitted nuisance predictions from \code{didml_nuisance()}
#'   with `estimand = "tc"` or `"both"`.
#'
#' @return List with `estimate`, `variance`, `scores`, `num_i`, `den_i`.
#'   Returns NA estimate with warning if first stage is too weak.
#'
#' @export
didml_tc <- function(W, pred) {
  N <- nrow(W)
  Y <- W$Y; D <- W$D; G <- W$G; Ti <- W$Ti
  GT <- G * Ti
  p_11 <- mean(GT)

  if (p_11 < 1e-10) {
    warning("No observations in treated group at post-period (G=1, T=1).", call. = FALSE)
    return(list(estimate = NA_real_, variance = NA_real_, scores = rep(NA_real_, N),
                num_i = rep(NA_real_, N), den_i = rep(NA_real_, N)))
  }

  m_Y_10 <- pred$m_Y_10; m_D_10 <- pred$m_D_10

  # Time trends by treatment status in control group
  delta_1 <- pred$mu_Y_101 - pred$mu_Y_100
  delta_0 <- pred$mu_Y_001 - pred$mu_Y_000

  # TC moment components
  g_Y <- Y - m_Y_10 - delta_1 * m_D_10 - (1 - m_D_10) * delta_0
  g_D <- D - m_D_10

  # Propensities and weights
  pi_11 <- pred$pi_11; pi_10 <- pred$pi_10
  w_10 <- G * (1 - Ti) * pi_11 / (pi_10 * p_11)
  w_101 <- D * (1 - G) * Ti * pi_11 / (pmax(pred$pi_101, .MIN_PROPENSITY) * p_11)
  w_100 <- D * (1 - G) * (1 - Ti) * pi_11 / (pmax(pred$pi_100, .MIN_PROPENSITY) * p_11)
  w_001 <- (1 - D) * (1 - G) * Ti * pi_11 / (pmax(pred$pi_001, .MIN_PROPENSITY) * p_11)
  w_000 <- (1 - D) * (1 - G) * (1 - Ti) * pi_11 / (pmax(pred$pi_000, .MIN_PROPENSITY) * p_11)

  # Set non-finite weights to NA
  for (nm in c("w_10", "w_101", "w_100", "w_001", "w_000")) {
    w <- get(nm)
    if (any(!is.finite(w))) {
      warning("Non-finite weights in ", nm, ". Affected obs set to NA.", call. = FALSE)
      w[!is.finite(w)] <- NA_real_
      assign(nm, w)
    }
  }

  # Residuals
  r_Y_10 <- Y - m_Y_10; r_D_10 <- D - m_D_10
  r_101 <- Y - pred$mu_Y_101; r_100 <- Y - pred$mu_Y_100
  r_001 <- Y - pred$mu_Y_001; r_000 <- Y - pred$mu_Y_000

  # Closed-form GMM
  num_i <- (GT / p_11) * g_Y -
    w_10 * (r_Y_10 - r_D_10 * (delta_0 - delta_1)) -
    m_D_10 * (w_101 * r_101 - w_100 * r_100) -
    (1 - m_D_10) * (w_001 * r_001 - w_000 * r_000)

  den_i <- (GT / p_11) * g_D - w_10 * r_D_10

  sum_den <- sum(den_i, na.rm = TRUE)
  if (abs(sum_den) < 1e-10) {
    warning("Weak first stage: TC denominator near zero. Returning NA.", call. = FALSE)
    return(list(estimate = NA_real_, variance = NA_real_, scores = rep(NA_real_, N),
                num_i = num_i, den_i = den_i))
  }

  theta_hat <- sum(num_i, na.rm = TRUE) / sum_den

  # Scores at theta_hat
  psi_i <- (GT / p_11) * (g_Y - g_D * theta_hat) -
    w_10 * (r_Y_10 - r_D_10 * (delta_0 - delta_1 + theta_hat)) -
    m_D_10 * (w_101 * r_101 - w_100 * r_100) -
    (1 - m_D_10) * (w_001 * r_001 - w_000 * r_000)

  G_hat <- -mean(den_i, na.rm = TRUE)
  ok <- is.finite(psi_i)
  V_hat <- mean(psi_i[ok]^2) / (G_hat^2)

  list(
    estimate = theta_hat,
    variance = V_hat,
    scores = psi_i,
    num_i = num_i,
    den_i = den_i
  )
}


#' Compute DML-Chang Scores (Sharp DID)
#'
#' Computes the doubly robust score for the sharp DID ATT estimand following
#' Chang (2020), equation 3.2. Designed for repeated cross-section settings
#' without an instrumental variable / fuzzy compliance.
#'
#' @details
#' Chang's (2020) influence function for the ATT in a 2x2 DID with
#' repeated cross-sections. In our notation (Chang's D = our G,
#' Chang's T = our Ti):
#'
#' \deqn{
#'   w_i = \frac{G_i - \hat{e}(X_i)}{p_T (1-p_T) \cdot p_0 \cdot (1 - \hat{e}(X_i))}
#' }
#' \deqn{
#'   \psi_i = w_i \left[(T_i - p_T) Y_i - \hat{\ell}_{20}(X_i)\right]
#' }
#' \deqn{
#'   \widehat{ATT} = \frac{1}{N} \sum_{i=1}^N \psi_i
#' }
#'
#' where \eqn{p_0 = \Pr(G=0) = 1 - \bar{G}}, \eqn{p_T = \bar{T}_i},
#' and \eqn{\ell_{20}(X) = E[(T_i - p_T) Y_i | X_i, G_i = 0]}.
#'
#' @param W Data frame with columns Y, D, G, Ti.
#' @param pred List of cross-fitted nuisance predictions from
#'   \code{didml_nuisance(iv = FALSE)}, containing \code{pG_raw},
#'   \code{ell_20}, and \code{pT}.
#'
#' @return List with \code{estimate}, \code{variance}, \code{scores},
#'   \code{num_i}, \code{den_i}. \code{den_i} is \code{rep(1, N)} for
#'   compatibility with the inference pipeline (sharp ATT has no ratio).
#'
#' @references
#' Chang, N.-C. (2020). Double/debiased machine learning for
#' difference-in-differences treatment effects.
#' \emph{The Econometrics Journal}, 23(2), 177-191.
#'
#' @export
didml_chang <- function(W, pred) {
  N <- nrow(W)
  Y <- W$Y; G <- W$G; Ti <- W$Ti

  pG_hat <- pred$pG_raw
  ell_20 <- pred$ell_20
  pT <- pred$pT

  # Marginal probabilities

  pG_bar <- mean(G)  # Pr(G = 1) = P(D=1) in Chang's notation

  if (sum(G) < 10 || sum(1 - G) < 10) {
    warning("Fewer than 10 treated or control observations. Cannot estimate ATT.", call. = FALSE)
    return(list(estimate = NA_real_, variance = NA_real_, scores = rep(NA_real_, N),
                num_i = rep(NA_real_, N), den_i = rep(1, N)))
  }

  # Weight: w_i = (G - e(X)) / (pT * (1-pT) * pG_bar * (1 - e(X)))
  # Chang (2020) eq. 3.2: denominator normalizes by P(D=1) * (1 - g_0(X))
  denom_w <- pT * (1 - pT) * pG_bar * (1 - pG_hat)

  # Guard against division by near-zero
  denom_w <- pmax(denom_w, 1e-10)

  w_i <- (G - pG_hat) / denom_w

  # Guard non-finite weights
  if (any(!is.finite(w_i))) {
    warning("Non-finite Chang weights (near-zero propensity). ",
            "Affected obs set to NA.", call. = FALSE)
    w_i[!is.finite(w_i)] <- NA_real_
  }

  # Score: psi_i = w_i * ((Ti - pT)*Y - ell_20(X))
  psi_i <- w_i * ((Ti - pT) * Y - ell_20)

  # ATT estimate
  ok <- is.finite(psi_i)
  ATT <- mean(psi_i[ok])

  # Variance: V = E[(psi_i - ATT)^2]
  V_hat <- mean((psi_i[ok] - ATT)^2)

  # den_i = 1 for all obs (no ratio estimand in sharp DID)
  den_i <- rep(1, N)

  list(
    estimate = ATT,
    variance = V_hat,
    scores = psi_i,
    num_i = psi_i,
    den_i = den_i
  )
}
