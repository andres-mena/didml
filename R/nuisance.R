#' Cross-Fitted Nuisance Estimation
#'
#' Estimates all nuisance functions required for DML estimation using
#' K-fold cross-fitting. Supports Lasso, Random Forest, neural networks,
#' OLS, or a user-supplied function.
#'
#' @details
#' When \code{iv = TRUE} (fuzzy DID), estimates conditional expectations
#' of Y and D within each (g,t) cell, plus a cross-fitted propensity score.
#' Joint cell probabilities are computed as Pr(G=g|X) * Pr(T=t),
#' assuming time period T is independent of group and covariates.
#' This holds in repeated cross-section designs where T indexes calendar
#' time. It does NOT hold in balanced panel data.
#'
#' When \code{iv = FALSE} (sharp DID / Chang 2020), estimates:
#' \enumerate{
#'   \item Cross-fitted propensity score e(X) = Pr(G=1|X).
#'   \item The pseudo-outcome regression ell_20(X) = E\[(Ti - pT)*Y \| X, G=0\]
#'     via cross-fitting within the control group (G=0).
#' }
#' No cell-specific conditional expectations are needed.
#'
#' When a cell-fold intersection has fewer than 10 training observations,
#' the function returns NA for those predictions (never imputes cell means).
#' Downstream score functions will produce NA estimates for observations
#' with missing nuisance predictions.
#'
#' @param Y Numeric vector of outcomes.
#' @param D Numeric vector of treatment.
#' @param G Binary vector of group assignment (0/1).
#' @param Ti Binary vector of time period (0/1).
#' @param X Numeric matrix of covariates.
#' @param method Character string (`"lasso"`, `"rf"`, `"nn"`, `"ols"`) or a
#'   function `f(Y_train, X_train, X_predict)` returning predicted values.
#' @param K Integer >= 2, number of cross-fitting folds (default 5).
#' @param iv Logical. If \code{TRUE} (default), estimate nuisance for fuzzy
#'   DID (Wald/TC). If \code{FALSE}, estimate nuisance for sharp DID
#'   (Chang 2020).
#' @param estimand Which estimand to prepare nuisance for. When \code{iv = TRUE}:
#'   `"wald"`, `"tc"`, or `"both"` (default). When \code{iv = FALSE}:
#'   `"chang"` (only valid option, set automatically).
#' @param seed Random seed for fold assignment (default `NULL`).
#'
#' @return A list with cross-fitted predictions for every observation.
#'   Entries may contain NA where estimation failed (small cells, collinearity).
#'   When \code{iv = FALSE}, the list contains: \code{pG_raw}, \code{pi_11},
#'   \code{pi_10}, \code{pi_01}, \code{pi_00}, \code{ell_20}, \code{pT},
#'   \code{folds}, and \code{estimand = "chang"}.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' N <- 500
#' X <- matrix(rnorm(N * 5), ncol = 5)
#' G <- rbinom(N, 1, 0.5)
#' Ti <- rbinom(N, 1, 0.5)
#' D <- rbinom(N, 1, plogis(X[,1] + G))
#' Y <- 1 + X[,1] + D + rnorm(N)
#' # Fuzzy DID nuisance
#' nuis_fuzzy <- didml_nuisance(Y, D, G, Ti, X, method = "ols", K = 2, iv = TRUE)
#' # Sharp DID nuisance (Chang 2020)
#' nuis_sharp <- didml_nuisance(Y, D, G, Ti, X, method = "ols", K = 2, iv = FALSE)
#' }
#'
#' @export
didml_nuisance <- function(Y, D, G, Ti, X,
                            method = "lasso",
                            K = 5L,
                            iv = TRUE,
                            estimand = "both",
                            seed = NULL) {
  # Validate method
  if (is.character(method)) {
    method <- match.arg(method, c("lasso", "rf", "nn", "ols"))
  } else if (!is.function(method)) {
    stop("`method` must be a character string or a function.", call. = FALSE)
  }
  if (!is.numeric(K) || length(K) != 1 || K < 2L)
    stop("`K` must be a single integer >= 2.", call. = FALSE)
  K <- as.integer(K)

  X <- as.matrix(X)
  N <- length(Y)
  .validate_inputs(Y, D, G, Ti, X)

  folds <- .make_folds(N, K, seed = seed)
  pT <- mean(Ti)

  # ---- Sharp DID path (Chang 2020) ----
  if (!iv) {
    estimand <- "chang"

    # Cross-fitted propensity score P(G=1|X)
    pG_raw <- rep(NA_real_, N)
    for (k in seq_len(K)) {
      train_k <- which(folds != k)
      pred_k <- which(folds == k)
      pG_raw[pred_k] <- .fit_nuisance_cell(G[train_k], X[train_k, , drop = FALSE],
                                             X[pred_k, , drop = FALSE], method,
                                             family = "binomial")
    }
    # Pre-clip to avoid numerical issues in downstream weight computation.
    # The trimming step in didml() applies the user's chosen rule on top.
    pG_raw <- pmax(pmin(pG_raw, 1 - .MIN_PROPENSITY), .MIN_PROPENSITY)

    # Joint cell probabilities
    pi_11 <- pG_raw * pT
    pi_10 <- pG_raw * (1 - pT)
    pi_01 <- (1 - pG_raw) * pT
    pi_00 <- (1 - pG_raw) * (1 - pT)

    # Pseudo-outcome: (Ti - pT) * Y, regressed on X within G=0
    pseudo_Y <- (Ti - pT) * Y
    idx_g0 <- which(G == 0)

    ell_20 <- rep(NA_real_, N)
    for (k in seq_len(K)) {
      train_k <- idx_g0[folds[idx_g0] != k]
      pred_k <- which(folds == k)
      ell_20[pred_k] <- .fit_nuisance_cell(pseudo_Y[train_k],
                                             X[train_k, , drop = FALSE],
                                             X[pred_k, , drop = FALSE], method)
    }

    result <- list(
      pG_raw = pG_raw,
      pi_11 = pi_11, pi_10 = pi_10, pi_01 = pi_01, pi_00 = pi_00,
      ell_20 = ell_20,
      pT = pT,
      folds = folds,
      estimand = "chang"
    )

    # Warn about NAs
    core_preds <- c("pG_raw", "ell_20")
    na_counts <- vapply(result[core_preds], function(x) sum(is.na(x)), integer(1L))
    if (any(na_counts > 0)) {
      bad <- core_preds[na_counts > 0]
      warning("Nuisance predictions contain NA in: ",
              paste(bad, collapse = ", "),
              ". Check cell sizes and method convergence.", call. = FALSE)
    }

    return(result)
  }

  # ---- Fuzzy DID path (iv = TRUE) ----
  estimand <- match.arg(estimand, c("wald", "tc", "both"))

  # Cell indices
  cells <- list(
    "10" = which(G == 1 & Ti == 0),
    "01" = which(G == 0 & Ti == 1),
    "00" = which(G == 0 & Ti == 0)
  )

  # Initialize prediction vectors as NA — never impute
  m_Y <- m_D <- list()
  for (key in names(cells)) {
    m_Y[[key]] <- rep(NA_real_, N)
    m_D[[key]] <- rep(NA_real_, N)
  }

  # Cross-fitted conditional expectations within each (g,t) cell
  for (key in names(cells)) {
    idx <- cells[[key]]
    for (k in seq_len(K)) {
      train_k <- idx[folds[idx] != k]
      pred_k <- which(folds == k)

      # .fit_nuisance_cell returns NA when cell too small or estimation fails
      m_Y[[key]][pred_k] <- .fit_nuisance_cell(Y[train_k], X[train_k, , drop = FALSE],
                                                 X[pred_k, , drop = FALSE], method)
      m_D[[key]][pred_k] <- .fit_nuisance_cell(D[train_k], X[train_k, , drop = FALSE],
                                                 X[pred_k, , drop = FALSE], method)
    }
  }

  # Cross-fitted propensity score P(G=1|X)
  pG_raw <- rep(NA_real_, N)
  for (k in seq_len(K)) {
    train_k <- which(folds != k)
    pred_k <- which(folds == k)
    pG_raw[pred_k] <- .fit_nuisance_cell(G[train_k], X[train_k, , drop = FALSE],
                                           X[pred_k, , drop = FALSE], method,
                                           family = "binomial")
  }
  # Pre-clip to avoid numerical issues in downstream weight computation
  pG_raw <- pmax(pmin(pG_raw, 1 - .MIN_PROPENSITY), .MIN_PROPENSITY)

  # Joint cell probabilities: pi_{gt}(X) = Pr(G=g|X) * Pr(T=t)
  # Assumes T independent of (G, X) — valid in repeated cross-sections
  pi_11 <- pG_raw * pT
  pi_10 <- pG_raw * (1 - pT)
  pi_01 <- (1 - pG_raw) * pT
  pi_00 <- (1 - pG_raw) * (1 - pT)

  result <- list(
    m_Y_10 = m_Y[["10"]], m_Y_01 = m_Y[["01"]], m_Y_00 = m_Y[["00"]],
    m_D_10 = m_D[["10"]], m_D_01 = m_D[["01"]], m_D_00 = m_D[["00"]],
    pG_raw = pG_raw,
    pi_11 = pi_11, pi_10 = pi_10, pi_01 = pi_01, pi_00 = pi_00,
    folds = folds
  )

  # Warn about NAs in core predictions
  core_preds <- c("m_Y_10", "m_Y_01", "m_Y_00", "m_D_10", "m_D_01", "m_D_00", "pG_raw")
  na_counts <- vapply(result[core_preds], function(x) sum(is.na(x)), integer(1L))
  if (any(na_counts > 0)) {
    bad <- core_preds[na_counts > 0]
    warning("Nuisance predictions contain NA in: ",
            paste(bad, collapse = ", "),
            ". Check cell sizes and method convergence.", call. = FALSE)
  }

  # TC-specific: treatment-conditional expectations in control group
  if (estimand %in% c("tc", "both")) {
    tc_cells <- list(
      "101" = which(D == 1 & G == 0 & Ti == 1),
      "100" = which(D == 1 & G == 0 & Ti == 0),
      "001" = which(D == 0 & G == 0 & Ti == 1),
      "000" = which(D == 0 & G == 0 & Ti == 0)
    )

    for (key in names(tc_cells)) {
      idx <- tc_cells[[key]]
      mu_key <- paste0("mu_Y_", key)
      pi_key <- paste0("pi_", key)
      result[[mu_key]] <- rep(NA_real_, N)
      result[[pi_key]] <- rep(NA_real_, N)

      for (k in seq_len(K)) {
        train_k <- idx[folds[idx] != k]
        pred_k <- which(folds == k)

        # mu_Y_{dgt}(X_i) trained on (d,0,t) cell, predicted for ALL obs.
        # Cell-membership indicators in the score zero out irrelevant terms.
        result[[mu_key]][pred_k] <- .fit_nuisance_cell(
          Y[train_k], X[train_k, , drop = FALSE],
          X[pred_k, , drop = FALSE], method)

        # Propensity for this (d,g,t) cell
        L_ind <- as.numeric(seq_len(N) %in% idx)
        result[[pi_key]][pred_k] <- .fit_nuisance_cell(
          L_ind[folds != k], X[folds != k, , drop = FALSE],
          X[pred_k, , drop = FALSE], method, family = "binomial")
      }
      result[[pi_key]] <- pmax(result[[pi_key]], .MIN_PROPENSITY)
    }
  }

  result
}
