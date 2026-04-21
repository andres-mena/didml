# =============================================================================
# 03_estimators.R — Stage 3: Compute estimators from saved nuisance (CHEAP)
#
# Loads nuisance predictions + raw data from Stage 2, computes all 7 estimators.
# DR-Wald and DR-TC use per-replication data-driven trimming (Remark 2 of
# Theorem 2: one-sided, minimize R_hat(alpha) on a grid). Does NOT require
# knowing tau. Non-DR estimators use PG_TRIM_DEFAULT (symmetric).
#
# Input:  data/simulations/nuisance/nuisance_n{N}.rds
# Output: output/simulations/estimates/estimates_n{N}.rds
#
#   Structure per file:
#   $estimates    — M x (7*4) matrix: {est, se, cov, ci_len} per estimator
#   $alpha_chosen — M x 1 vector: data-driven alpha per replication
#   $trim_used    — list describing trimming method
#   $n, $M, $tau, $estimators, $timestamp
#
# Usage:
#   Rscript code/R/sim/03_estimators.R                       # data-driven (default)
#   Rscript code/R/sim/03_estimators.R --n 500               # single n
#   Rscript code/R/sim/03_estimators.R --trim 0.15           # override: fixed symmetric
#   Rscript code/R/sim/03_estimators.R --uniform-trim        # PG_TRIM_DEFAULT symmetric
#
# =============================================================================

source("code/R/sim/00_config.R")
source("code/R/04_estimators.R")  # compute_dd_wald, compute_dd_tc

args <- commandArgs(trailingOnly = TRUE)
n_filter <- NULL
trim_override <- NULL
ML_METHOD <- "lasso"
USE_DATA_DRIVEN <- TRUE  # default: per-replication data-driven trimming
for (i in seq_along(args)) {
  if (args[i] == "--n" && i < length(args)) n_filter <- as.integer(args[i + 1])
  if (args[i] == "--trim" && i < length(args)) {
    trim_override <- as.numeric(args[i + 1])
    USE_DATA_DRIVEN <- FALSE
  }
  if (args[i] == "--uniform-trim") USE_DATA_DRIVEN <- FALSE
  if (args[i] == "--method" && i < length(args)) ML_METHOD <- args[i + 1]
}
n_run <- if (!is.null(n_filter)) n_filter else N_VALUES

# =============================================================================
# DATA-DRIVEN TRIMMING: Remark 2 of Theorem 2 (Corollary 2.1, one-sided)
#
# Minimize R_hat(alpha) over alpha in {0, 0.01, ..., 0.50}.
#
# Full Corollary 2.1 formula (general g):
#   R_hat(alpha) = sum[e/(1-e) * 1{e<=1-a}] / (sum[e * g_hat(e) * 1{e<=1-a}])^2
#
# where g_hat(e) = E[DID_D | e(X)=e] is estimated via local linear regression
# of the cross-fitted DID_D on the estimated propensity score.
#
# Under first-stage independence (Cor 2.2): g_hat = constant, denominator
# simplifies to (sum[e * 1])^2. The general formula nests this as a special case.
#
# Does NOT use tau. Uses cross-fitted pG_raw and m_D_{gt} from Stage 2.
# =============================================================================
ALPHA_GRID_DD <- seq(0, 0.50, by = 0.01)

compute_data_driven_alpha <- function(pG_raw, DID_D_hat) {
  n <- length(pG_raw)

  # Estimate g_hat(e) = E[DID_D | e(X) = e] via local linear regression
  # Use loess with span adapted to sample size
  span_val <- max(0.3, min(0.75, 50 / n))
  g_fit <- tryCatch({
    loess(DID_D_hat ~ pG_raw, span = span_val, degree = 1)
  }, error = function(e) NULL)

  if (!is.null(g_fit)) {
    g_hat <- pmax(predict(g_fit, newdata = data.frame(pG_raw = pG_raw)), 0.001)
    g_hat[is.na(g_hat)] <- mean(DID_D_hat, na.rm = TRUE)
  } else {
    # Fallback: constant g (Corollary 2.2)
    g_hat <- rep(mean(DID_D_hat), n)
  }

  best_alpha <- 0
  best_R <- Inf

  for (alpha in ALPHA_GRID_DD) {
    keep <- pG_raw <= (1 - alpha)
    if (sum(keep) < 10) next

    numerator   <- sum(pG_raw[keep] / (1 - pG_raw[keep]))
    denominator <- sum(pG_raw[keep] * g_hat[keep])^2
    if (denominator < 1e-10) next

    R_val <- numerator / denominator
    if (R_val < best_R) {
      best_R <- R_val
      best_alpha <- alpha
    }
  }
  best_alpha
}

# =============================================================================
# COMPUTE ALL 7 ESTIMATORS FROM SAVED NUISANCE
# =============================================================================
compute_all_from_nuisance <- function(nuis, tau, alpha_dd, fixed_trim, use_dd) {
  Y <- nuis$Y; D <- nuis$D; G <- nuis$G; Ti <- nuis$Ti
  n <- length(Y)
  GT <- G * Ti; p_11 <- mean(GT); pT <- mean(Ti)

  if (p_11 < 0.01) return(rep(NA_real_, 14))

  out <- numeric(14)
  names(out) <- paste0(rep(ESTIMATOR_NAMES, each = 2), c("_est", "_se"))

  # --- Helper: Wald ratio with SE ---
  # Non-DR estimators always use symmetric clipping at fixed_trim
  wald_ratio <- function(mY, mD, pG_raw, use_trim) {
    DID_Y <- Y - mY[["10"]] - mY[["01"]] + mY[["00"]]
    DID_D <- D - mD[["10"]] - mD[["01"]] + mD[["00"]]
    den <- sum(GT * DID_D)
    if (abs(den) < 1e-10) return(c(NA, NA))
    theta <- sum(GT * DID_Y) / den

    pG <- pmax(pmin(pG_raw, 1 - use_trim), use_trim)
    pi11 <- pG * pT; pi10 <- pG * (1 - pT)
    pi01 <- (1 - pG) * pT; pi00 <- (1 - pG) * (1 - pT)
    w10 <- G * (1 - Ti) * pi11 / (pi10 * p_11)
    w01 <- (1 - G) * Ti * pi11 / (pi01 * p_11)
    w00 <- (1 - G) * (1 - Ti) * pi11 / (pi00 * p_11)

    eps <- Y - theta * D
    me10 <- mY[["10"]] - theta * mD[["10"]]
    me01 <- mY[["01"]] - theta * mD[["01"]]
    me00 <- mY[["00"]] - theta * mD[["00"]]
    psi <- (GT / p_11) * (eps - me10 - me01 + me00) -
      w10 * (eps - me10) - w01 * (eps - me01) + w00 * (eps - me00)
    Gh <- -mean((GT / p_11) * DID_D -
                  w10 * (D - mD[["10"]]) - w01 * (D - mD[["01"]]) + w00 * (D - mD[["00"]]))
    if (abs(Gh) < 1e-10) return(c(NA, NA))
    c(theta, sqrt(mean(psi^2) / (Gh^2) / n))
  }

  # --- Helper: TC ratio with SE (full CD'H Eq 69 IF for variance) ---
  # Point estimate uses g^tc only (plug-in, no FSIF).
  # SE uses the full efficient IF from CD'H supplement Eq 69,
  # which includes FSIF correction terms for all cells.
  # Without these terms, the SE underestimates the true SD by 3-4x.
  tc_ratio <- function(mY10, mD10, mu_Y, pi_tc, pG_raw, use_trim) {
    d1 <- mu_Y[["101"]] - mu_Y[["100"]]
    d0 <- mu_Y[["001"]] - mu_Y[["000"]]
    gY <- Y - mY10 - d1 * mD10 - (1 - mD10) * d0
    gD <- D - mD10
    den <- sum(GT * gD)
    if (abs(den) < 1e-10) return(c(NA, NA))
    theta <- sum(GT * gY) / den

    # Full TC IF for variance estimation (CD'H Eq 69)
    pG <- pmax(pmin(pG_raw, 1 - use_trim), use_trim)
    pi11 <- pG * pT; pi10 <- pG * (1 - pT)

    # w_10 correction (same structure as Wald)
    w_10 <- G * (1 - Ti) * pi11 / (pi10 * p_11)
    r_Y_10 <- Y - mY10; r_D_10 <- D - mD10

    # TC cell corrections
    pi_101 <- pmax(pi_tc[["101"]], 0.01)
    pi_100 <- pmax(pi_tc[["100"]], 0.01)
    pi_001 <- pmax(pi_tc[["001"]], 0.01)
    pi_000 <- pmax(pi_tc[["000"]], 0.01)

    w_101 <- D * (1 - G) * Ti * pi11 / (pi_101 * p_11)
    w_100 <- D * (1 - G) * (1 - Ti) * pi11 / (pi_100 * p_11)
    w_001 <- (1 - D) * (1 - G) * Ti * pi11 / (pi_001 * p_11)
    w_000 <- (1 - D) * (1 - G) * (1 - Ti) * pi11 / (pi_000 * p_11)

    r_101 <- Y - mu_Y[["101"]]; r_100 <- Y - mu_Y[["100"]]
    r_001 <- Y - mu_Y[["001"]]; r_000 <- Y - mu_Y[["000"]]

    # Full augmented score psi^tc = g^tc + phi^tc
    psi <- (GT / p_11) * (gY - gD * theta) -
      w_10 * (r_Y_10 - r_D_10 * (d0 - d1 + theta)) -
      mD10 * (w_101 * r_101 - w_100 * r_100) -
      (1 - mD10) * (w_001 * r_001 - w_000 * r_000)

    Gh <- -mean((GT / p_11) * gD - w_10 * r_D_10)
    if (abs(Gh) < 1e-10) return(c(NA, NA))
    c(theta, sqrt(mean(psi^2) / (Gh^2) / n))
  }

  # ---- 1. Wald (cell means, no covariates) ----
  mY_cell <- mD_cell <- list()
  for (g in 0:1) for (t in 0:1) {
    k <- paste0(g, t); idx <- G == g & Ti == t
    mY_cell[[k]] <- rep(mean(Y[idx]), n)
    mD_cell[[k]] <- rep(mean(D[idx]), n)
  }
  out[1:2] <- wald_ratio(mY_cell, mD_cell, rep(mean(G), n), 0.01)

  # ---- 2. Wald_X Lasso ----
  mY_l <- list("10" = nuis$fs_lasso_Y_10, "01" = nuis$fs_lasso_Y_01, "00" = nuis$fs_lasso_Y_00)
  mD_l <- list("10" = nuis$fs_lasso_D_10, "01" = nuis$fs_lasso_D_01, "00" = nuis$fs_lasso_D_00)
  out[3:4] <- wald_ratio(mY_l, mD_l, nuis$fs_pG, fixed_trim)

  # ---- 3. Wald_X Sieve ----
  mY_s <- list("10" = nuis$fs_sieve_Y_10, "01" = nuis$fs_sieve_Y_01, "00" = nuis$fs_sieve_Y_00)
  mD_s <- list("10" = nuis$fs_sieve_D_10, "01" = nuis$fs_sieve_D_01, "00" = nuis$fs_sieve_D_00)
  out[5:6] <- wald_ratio(mY_s, mD_s, nuis$fs_pG, fixed_trim)

  # ---- 4. TC_X Lasso ----
  mu_l <- list("101" = nuis$fs_lasso_mu_101, "100" = nuis$fs_lasso_mu_100,
               "001" = nuis$fs_lasso_mu_001, "000" = nuis$fs_lasso_mu_000)
  pi_tc_l <- list("101" = nuis$fs_pi_101, "100" = nuis$fs_pi_100,
                   "001" = nuis$fs_pi_001, "000" = nuis$fs_pi_000)
  out[7:8] <- tc_ratio(nuis$fs_lasso_Y_10, nuis$fs_lasso_D_10, mu_l,
                        pi_tc_l, nuis$fs_pG, fixed_trim)

  # ---- 5. TC_X Sieve ----
  mu_s <- list("101" = nuis$fs_sieve_mu_101, "100" = nuis$fs_sieve_mu_100,
               "001" = nuis$fs_sieve_mu_001, "000" = nuis$fs_sieve_mu_000)
  out[9:10] <- tc_ratio(nuis$fs_sieve_Y_10, nuis$fs_sieve_D_10, mu_s,
                         pi_tc_l, nuis$fs_pG, fixed_trim)

  # ---- 6. DR-Wald (cross-fitted + FSIF) ----
  W_df <- data.frame(Y = Y, D = D, G = G, Ti = Ti)

  # One-sided trimming: pG <- min(pG, 1 - alpha), floor at 0.001
  cf_w <- nuis$cf
  if (use_dd) {
    pG_w <- pmin(cf_w$pG_raw, 1 - alpha_dd)
    pG_w <- pmax(pG_w, 0.001)
  } else {
    pG_w <- pmax(pmin(cf_w$pG_raw, 1 - fixed_trim), fixed_trim)
  }
  cf_w$pi_11 <- pG_w * pT
  cf_w$pi_10 <- pG_w * (1 - pT)
  cf_w$pi_01 <- (1 - pG_w) * pT
  cf_w$pi_00 <- (1 - pG_w) * (1 - pT)

  drw <- tryCatch(compute_dd_wald(W_df, cf_w), error = function(e) NULL)
  out[11:12] <- if (!is.null(drw)) c(drw$estimate, drw$se) else c(NA, NA)

  # ---- 7. DR-TC (cross-fitted + FSIF, same data-driven alpha) ----
  cf_t <- nuis$cf
  if (use_dd) {
    pG_t <- pmin(cf_t$pG_raw, 1 - alpha_dd)
    pG_t <- pmax(pG_t, 0.001)
  } else {
    pG_t <- pmax(pmin(cf_t$pG_raw, 1 - fixed_trim), fixed_trim)
  }
  cf_t$pi_11 <- pG_t * pT
  cf_t$pi_10 <- pG_t * (1 - pT)
  cf_t$pi_01 <- (1 - pG_t) * pT
  cf_t$pi_00 <- (1 - pG_t) * (1 - pT)

  drt <- tryCatch(compute_dd_tc(W_df, cf_t), error = function(e) NULL)
  out[13:14] <- if (!is.null(drt)) c(drt$estimate, drt$se) else c(NA, NA)

  out
}

# =============================================================================
# MAIN LOOP
# =============================================================================
cat("=============================================================================\n")
cat("03_estimators.R — Compute estimators from saved nuisance\n")
if (USE_DATA_DRIVEN) {
  cat("  Trimming: DATA-DRIVEN per replication (Remark 2, one-sided)\n")
  cat(sprintf("  Alpha grid: %.2f to %.2f (step %.2f)\n",
              min(ALPHA_GRID_DD), max(ALPHA_GRID_DD), diff(ALPHA_GRID_DD)[1]))
} else if (!is.null(trim_override)) {
  cat(sprintf("  Trimming: FIXED SYMMETRIC at %.2f (--trim override)\n", trim_override))
} else {
  cat(sprintf("  Trimming: UNIFORM SYMMETRIC at %.2f (default)\n", PG_TRIM_DEFAULT))
}
cat(sprintf("  n: %s\n", paste(n_run, collapse = ",")))
cat("=============================================================================\n\n")

for (n_val in n_run) {
  infile <- nuisance_file(n_val, ML_METHOD)
  if (!file.exists(infile)) {
    cat(sprintf("SKIP n=%d: %s not found. Run 02_nuisance.R first.\n", n_val, infile))
    next
  }

  # Fixed trim for non-DR estimators (and fallback)
  fixed_trim <- if (!is.null(trim_override)) trim_override else PG_TRIM_DEFAULT

  cat(sprintf("[%s] n=%d: loading nuisance...\n", format(Sys.time(), "%H:%M:%S"), n_val))
  saved <- readRDS(infile)
  all_nuis <- saved$nuisance
  M_actual <- length(all_nuis)

  cat(sprintf("  Computing %d reps x %d estimators", M_actual, length(ESTIMATOR_NAMES)))
  if (USE_DATA_DRIVEN) cat(" (data-driven trim)")
  else cat(sprintf(" (fixed trim=%.2f)", fixed_trim))
  cat("...\n")
  t0 <- Sys.time()

  # Result matrix
  col_names <- paste0(rep(ESTIMATOR_NAMES, each = 4), c("_est", "_se", "_cov", "_ci_len"))
  results <- matrix(NA_real_, M_actual, length(col_names))
  colnames(results) <- col_names

  # Track data-driven alpha per replication
  alpha_vec <- rep(NA_real_, M_actual)

  for (r in 1:M_actual) {
    nuis <- all_nuis[[r]]
    if (inherits(nuis, "error") || !is.list(nuis) || is.null(nuis$Y)) next

    # Compute data-driven alpha for this replication
    # DID_D_hat = D - m_D_10(X) - m_D_01(X) + m_D_00(X) using cross-fitted predictions
    alpha_r <- if (USE_DATA_DRIVEN) {
      tryCatch({
        DID_D_hat <- nuis$D - nuis$cf$m_D_10 - nuis$cf$m_D_01 + nuis$cf$m_D_00
        compute_data_driven_alpha(nuis$cf$pG_raw, DID_D_hat)
      }, error = function(e) 0)
    } else 0
    alpha_vec[r] <- alpha_r

    raw <- tryCatch(
      compute_all_from_nuisance(nuis, TAU, alpha_r, fixed_trim, USE_DATA_DRIVEN),
      error = function(e) rep(NA_real_, 14))

    for (j in seq_along(ESTIMATOR_NAMES)) {
      nm <- ESTIMATOR_NAMES[j]
      e_val <- raw[2 * j - 1]; s_val <- raw[2 * j]
      results[r, paste0(nm, "_est")] <- e_val
      results[r, paste0(nm, "_se")]  <- s_val
      if (is.finite(e_val) && is.finite(s_val) && s_val > 0) {
        results[r, paste0(nm, "_cov")]    <- as.numeric((e_val - 1.96 * s_val) <= TAU &
                                                           TAU <= (e_val + 1.96 * s_val))
        results[r, paste0(nm, "_ci_len")] <- 2 * 1.96 * s_val
      }
    }
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  # Save
  outfile <- estimates_file(n_val, ML_METHOD)
  saveRDS(list(
    estimates     = results,
    alpha_chosen  = alpha_vec,
    trim_used     = list(
      mode       = if (USE_DATA_DRIVEN) "data-driven" else "fixed",
      method     = if (USE_DATA_DRIVEN) "Remark 2, one-sided, R_hat(alpha) grid" else "symmetric clipping",
      fixed_trim = fixed_trim,
      non_dr     = PG_TRIM_DEFAULT,
      tc_cell    = TC_TRIM
    ),
    n = n_val, M = M_actual, tau = TAU,
    estimators = ESTIMATOR_NAMES,
    timestamp  = Sys.time()
  ), outfile)

  mb <- round(file.size(outfile) / 1e6, 3)
  cat(sprintf("  Done in %.1fs → %s (%.3f MB)\n", elapsed, outfile, mb))

  # Print alpha summary
  if (USE_DATA_DRIVEN) {
    ok_a <- is.finite(alpha_vec)
    cat(sprintf("  Data-driven alpha: mean=%.3f, median=%.3f, sd=%.3f, range=[%.2f, %.2f]\n",
                mean(alpha_vec[ok_a]), median(alpha_vec[ok_a]), sd(alpha_vec[ok_a]),
                min(alpha_vec[ok_a]), max(alpha_vec[ok_a])))
  }

  # Print estimator summary
  for (nm in ESTIMATOR_NAMES) {
    e <- results[, paste0(nm, "_est")]
    ok <- is.finite(e)
    if (sum(ok) > 5) {
      cat(sprintf("    %15s: bias=%7.3f  sd=%6.3f  cov=%.2f  (%d valid)\n",
                  nm, mean(e[ok]) - TAU, sd(e[ok]),
                  mean(results[ok, paste0(nm, "_cov")], na.rm = TRUE), sum(ok)))
    }
  }
  cat("\n")
}

cat("=== Stage 3 complete ===\n")
