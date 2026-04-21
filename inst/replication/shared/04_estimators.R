# =============================================================================
# 04_estimators.R
# Purpose:  DD-Wald and DD-TC estimators (doubly robust GMM)
# Author:   Andres Mena (Brown University)
# Input:    Data W, cross-fitted predictions
# Output:   Point estimate, SE, CI
# =============================================================================

compute_dd_wald <- function(W, pred) {
  # DML-Wald estimator: locally robust GMM for the Wald-DID estimand
  #
  # The augmented score is:
  #   psi^wald = g^wald + phi^wald
  # where g^wald is the identifying moment and phi^wald is the FSIF correction.
  #
  # Sign convention: (-, -, +) for correction terms w_10, w_01, w_00.
  # This matches the DID structure DID = Y - m_10 - m_01 + m_00
  # and is equivalent to CD'H (2018) supplement p.45, ψ*_DID(g),
  # which has cell signs (+, -, -, +) for (1,1), (1,0), (0,1), (0,0).
  # Verified empirically: r = 1.0000 with CD'H p.45 IF.
  # See explorations/if_equivalence_report/ for full verification.
  #
  # Returns: list(estimate, se, ci_lower, ci_upper, scores)

  N <- nrow(W)
  Y <- W$Y
  D <- W$D
  G <- W$G
  Ti <- W$Ti
  GT <- G * Ti

  # Sample proportion p_11
  p_11 <- mean(GT)

  # Cross-fitted conditional expectations
  m_Y_10 <- pred$m_Y_10
  m_Y_01 <- pred$m_Y_01
  m_Y_00 <- pred$m_Y_00
  m_D_10 <- pred$m_D_10
  m_D_01 <- pred$m_D_01
  m_D_00 <- pred$m_D_00

  # DID of conditional expectations
  DID_Y <- Y - m_Y_10 - m_Y_01 + m_Y_00
  DID_D <- D - m_D_10 - m_D_01 + m_D_00

  # Naive Wald-DID (no debiasing, for comparison)
  wald_naive <- sum(GT * DID_Y) / sum(GT * DID_D)

  # Cross-fitted propensity scores
  pi_10 <- pred$pi_10
  pi_01 <- pred$pi_01
  pi_00 <- pred$pi_00
  pi_11 <- pred$pi_11

  # Riesz representer weights following CD'H (2018) Eq 68:
  # alpha_{gt}(X) = I(G=g,T=t) * E(GT|X) / (E(I_{gt}|X) * p_11)
  #               = I(G=g,T=t) * pi_11(X) / (pi_gt(X) * p_11)
  # Signs follow CD'H Eq 68: negative for (1,0), positive for (0,1), negative for (0,0)
  w_10 <- (G * (1 - Ti)) * pi_11 / (pi_10 * p_11)
  w_01 <- ((1 - G) * Ti) * pi_11 / (pi_01 * p_11)
  w_00 <- ((1 - G) * (1 - Ti)) * pi_11 / (pi_00 * p_11)

  # Residuals
  r_Y_10 <- Y - m_Y_10
  r_Y_01 <- Y - m_Y_01
  r_Y_00 <- Y - m_Y_00
  r_D_10 <- D - m_D_10
  r_D_01 <- D - m_D_01
  r_D_00 <- D - m_D_00

  # --- GMM estimation via closed-form ---
  # psi^wald = (GT/p_11)(DID_Y - DID_D*theta)
  #   - w_10*(r_Y_10 - r_D_10*theta)   [negative: m_10 has sign -1 in DID]
  #   - w_01*(r_Y_01 - r_D_01*theta)   [negative: m_01 has sign -1 in DID]
  #   + w_00*(r_Y_00 - r_D_00*theta)   [positive: m_00 has sign +1 in DID]
  #
  # FSIF correction signs match the signs of each m_gt in the DID:
  #   DID = Y - m_10 - m_01 + m_00  →  signs: -, -, +
  # This ensures double robustness: E[psi] = 0 for ANY m_hat when pi is correct.
  #
  # Setting E[psi] = 0: numerator / denominator = theta

  numerator_i <- (GT / p_11) * DID_Y -
    w_10 * r_Y_10 - w_01 * r_Y_01 + w_00 * r_Y_00

  denominator_i <- (GT / p_11) * DID_D -
    w_10 * r_D_10 - w_01 * r_D_01 + w_00 * r_D_00

  theta_hat <- sum(numerator_i) / sum(denominator_i)

  # --- Compute individual scores at theta_hat (for variance estimation) ---
  psi_i <- (GT / p_11) * (DID_Y - DID_D * theta_hat) -
    w_10 * (r_Y_10 - r_D_10 * theta_hat) -
    w_01 * (r_Y_01 - r_D_01 * theta_hat) +
    w_00 * (r_Y_00 - r_D_00 * theta_hat)

  # --- Variance estimation (Theorem 1) ---
  # V^wald = [G^wald]^{-2} * E[psi^2]
  # G^wald = -mean(denominator_i)
  G_hat <- -mean(denominator_i)

  V_hat <- mean(psi_i^2) / (G_hat^2)
  se_hat <- sqrt(V_hat / N)

  # 95% CI
  ci_lower <- theta_hat - 1.96 * se_hat
  ci_upper <- theta_hat + 1.96 * se_hat

  return(list(
    estimate = theta_hat,
    se = se_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    variance = V_hat,
    scores = psi_i,
    wald_naive = wald_naive
  ))
}


compute_wald_did_naive <- function(W, X_names = NULL) {
  # Naive Wald-DID with Lasso first step but NO debiasing (no FSIF correction)
  # This shows the bias introduced by ML regularization without orthogonalization.
  # Now also computes SE via the non-orthogonal score variance.

  N <- nrow(W)
  Y <- W$Y
  D <- W$D
  G <- W$G
  Ti <- W$Ti
  GT <- G * Ti

  if (is.null(X_names)) X_names <- grep("^X", names(W), value = TRUE)

  # Lasso within each (g,t) cell — same ML method, but NO cross-fitting, NO FSIF
  m_Y <- m_D <- list()
  for (g_val in 0:1) {
    for (t_val in 0:1) {
      key <- paste0(g_val, t_val)
      idx <- which(W$G == g_val & W$Ti == t_val)

      if (length(idx) < 10) {
        m_Y[[key]] <- rep(mean(W$Y[idx]), N)
        m_D[[key]] <- rep(mean(W$D[idx]), N)
        next
      }

      X_sub <- as.matrix(W[idx, X_names, drop = FALSE])
      X_all <- as.matrix(W[, X_names, drop = FALSE])

      # Lasso for Y
      fit_Y <- tryCatch({
        mod <- glmnet::cv.glmnet(X_sub, W$Y[idx], alpha = 1, nfolds = 5)
        as.vector(predict(mod, newx = X_all, s = "lambda.min"))
      }, error = function(e) rep(mean(W$Y[idx]), N))

      # Lasso for D
      fit_D <- tryCatch({
        mod <- glmnet::cv.glmnet(X_sub, W$D[idx], alpha = 1, nfolds = 5)
        as.vector(predict(mod, newx = X_all, s = "lambda.min"))
      }, error = function(e) rep(mean(W$D[idx]), N))

      m_Y[[key]] <- fit_Y
      m_D[[key]] <- fit_D
    }
  }

  DID_Y <- Y - m_Y[["10"]] - m_Y[["01"]] + m_Y[["00"]]
  DID_D <- D - m_D[["10"]] - m_D[["01"]] + m_D[["00"]]

  p_11 <- mean(GT)
  theta_hat <- sum(GT * DID_Y) / sum(GT * DID_D)

  # SE via CD'H efficient influence function (Eq 68 / Eq A29 from supplement)
  # The naive estimator uses g^wald for the POINT ESTIMATE (no orthogonalization),
  # but the CD'H IF for the VARIANCE (which is the correct asymptotic variance).
  # This requires propensity scores pi_{gt}(X) — estimated via logit-lasso on full sample.

  # Estimate propensity scores (no cross-fitting — matches the naive approach)
  pi_est <- list()
  X_all <- as.matrix(W[, X_names, drop = FALSE])
  for (g_val in 0:1) {
    for (t_val in 0:1) {
      key <- paste0(g_val, t_val)
      L_ind <- as.numeric(W$G == g_val & W$Ti == t_val)
      pi_est[[key]] <- tryCatch({
        mod <- glmnet::cv.glmnet(X_all, L_ind, alpha = 1, family = "binomial", nfolds = 5)
        pmax(pmin(as.vector(predict(mod, newx = X_all, s = "lambda.min", type = "response")), 0.90), 0.10)
      }, error = function(e) rep(mean(L_ind), N))
    }
  }

  # CD'H Eq 68: psi^X_DID (efficient IF with covariates)
  # eps = Y - theta*D;  m_eps_gt = m_Y_gt - theta*m_D_gt
  eps <- Y - theta_hat * D
  m_eps_10 <- m_Y[["10"]] - theta_hat * m_D[["10"]]
  m_eps_01 <- m_Y[["01"]] - theta_hat * m_D[["01"]]
  m_eps_00 <- m_Y[["00"]] - theta_hat * m_D[["00"]]

  # Riesz weights: I(G=g,T=t)*pi_11/(pi_gt*p_11)
  w_10 <- G * (1 - Ti) * pi_est[["11"]] / (pi_est[["10"]] * p_11)
  w_01 <- (1 - G) * Ti * pi_est[["11"]] / (pi_est[["01"]] * p_11)
  w_00 <- (1 - G) * (1 - Ti) * pi_est[["11"]] / (pi_est[["00"]] * p_11)

  # CD'H IF score at theta_hat
  DID_eps <- eps - m_eps_10 - m_eps_01 + m_eps_00
  psi_cdh_i <- (GT / p_11) * DID_eps -
    w_10 * (eps - m_eps_10) - w_01 * (eps - m_eps_01) + w_00 * (eps - m_eps_00)

  G_hat <- -mean((GT / p_11) * DID_D -
    w_10 * (D - m_D[["10"]]) - w_01 * (D - m_D[["01"]]) + w_00 * (D - m_D[["00"]]))
  V_hat <- mean(psi_cdh_i^2) / (G_hat^2)
  se_hat <- sqrt(V_hat / N)

  ci_lower <- theta_hat - 1.96 * se_hat
  ci_upper <- theta_hat + 1.96 * se_hat

  return(list(
    estimate = theta_hat,
    se = se_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}


# =============================================================================
# DML-TC: Locally Robust Time-Corrected Estimator
# =============================================================================

compute_dd_tc <- function(W, pred) {
  # DML-TC estimator: locally robust GMM for the Time-Corrected estimand
  # Equivalent to CD'H (2018) supplement p.45, ψ*_TC(g), up to normalization.
  # Verified empirically: r = 1.0000 with CD'H TC IF.
  # See explorations/if_equivalence_report/ for full verification.
  # Uses CD'H (2018) Eq 69 formulation with correct Riesz weights
  # Requires: cross-fitted mu^Y_{dgt} and pi_{dgt} in pred

  N <- nrow(W)
  Y <- W$Y; D <- W$D; G <- W$G; Ti <- W$Ti; GT <- G * Ti
  p_11 <- mean(GT)

  # Marginal conditional expectations (same as Wald)
  m_Y_10 <- pred$m_Y_10
  m_D_10 <- pred$m_D_10

  # Treatment-conditional expectations (TC-specific)
  mu_Y_101 <- pred$mu_Y_101  # E[Y|D=1,G=0,T=1,X]
  mu_Y_100 <- pred$mu_Y_100  # E[Y|D=1,G=0,T=0,X]
  mu_Y_001 <- pred$mu_Y_001  # E[Y|D=0,G=0,T=1,X]
  mu_Y_000 <- pred$mu_Y_000  # E[Y|D=0,G=0,T=0,X]

  # Conditional time trends by treatment status
  delta_1 <- mu_Y_101 - mu_Y_100  # trend for D=1 in control group
  delta_0 <- mu_Y_001 - mu_Y_000  # trend for D=0 in control group

  # Propensity scores
  pi_11 <- pred$pi_11; pi_10 <- pred$pi_10
  pi_101 <- pred$pi_101; pi_100 <- pred$pi_100
  pi_001 <- pred$pi_001; pi_000 <- pred$pi_000

  # --- TC moment function g^tc (Eq 6) ---
  # g^tc = (GT/p_11)[Y - m_Y_10 - delta_1*m_D_10 - (1-m_D_10)*delta_0 - (D-m_D_10)*theta]
  g_Y <- Y - m_Y_10 - delta_1 * m_D_10 - (1 - m_D_10) * delta_0
  g_D <- D - m_D_10

  # --- TC FSIF phi^tc (Eq 9, corrected weights) ---
  # Weight factor: pi_11/(p_11 * pi_{cell})
  # Term 1: -(G(1-T)*pi_11/(pi_10*p_11)) * (Y - m_Y_10 - (D-m_D_10)*(delta_1-delta_0+theta))
  w_10 <- G * (1 - Ti) * pi_11 / (pi_10 * p_11)

  # Terms 2-3: -m_D_10 * pi_11/p_11 * (D(1-G)T/pi_101*(Y-mu_101) - D(1-G)(1-T)/pi_100*(Y-mu_100))
  w_101 <- D * (1 - G) * Ti * pi_11 / (pi_101 * p_11)
  w_100 <- D * (1 - G) * (1 - Ti) * pi_11 / (pi_100 * p_11)

  # Terms 4-5: -(1-m_D_10) * pi_11/p_11 * ((1-D)(1-G)T/pi_001*(Y-mu_001) - (1-D)(1-G)(1-T)/pi_000*(Y-mu_000))
  w_001 <- (1 - D) * (1 - G) * Ti * pi_11 / (pi_001 * p_11)
  w_000 <- (1 - D) * (1 - G) * (1 - Ti) * pi_11 / (pi_000 * p_11)

  # Residuals for TC FSIF
  r_Y_10 <- Y - m_Y_10
  r_D_10 <- D - m_D_10
  r_101 <- Y - mu_Y_101
  r_100 <- Y - mu_Y_100
  r_001 <- Y - mu_Y_001
  r_000 <- Y - mu_Y_000

  # --- Closed-form: numerator and denominator ---
  # psi^tc = g^tc(theta) + phi^tc(theta)
  # g^tc numerator: (GT/p_11) * g_Y
  # phi^tc numerator (Y part): correction terms for Y residuals

  # For the linear solve, separate Y and D*theta components:
  # Numerator: sum of Y-related terms
  num_i <- (GT / p_11) * g_Y -
    w_10 * (r_Y_10 - r_D_10 * (delta_0 - delta_1)) -
    m_D_10 * (w_101 * r_101 - w_100 * r_100) -
    (1 - m_D_10) * (w_001 * r_001 - w_000 * r_000)

  # Denominator: sum of D-related terms (coefficient of theta)
  den_i <- (GT / p_11) * g_D -
    w_10 * r_D_10

  theta_hat <- sum(num_i) / sum(den_i)

  # --- Individual scores at theta_hat ---
  psi_i <- (GT / p_11) * (g_Y - g_D * theta_hat) -
    w_10 * (r_Y_10 - r_D_10 * (delta_0 - delta_1 + theta_hat)) -
    m_D_10 * (w_101 * r_101 - w_100 * r_100) -
    (1 - m_D_10) * (w_001 * r_001 - w_000 * r_000)

  # --- Variance estimation ---
  G_hat <- -mean(den_i)
  V_hat <- mean(psi_i^2) / (G_hat^2)
  se_hat <- sqrt(V_hat / N)

  ci_lower <- theta_hat - 1.96 * se_hat
  ci_upper <- theta_hat + 1.96 * se_hat

  return(list(
    estimate = theta_hat,
    se = se_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    variance = V_hat,
    scores = psi_i
  ))
}
