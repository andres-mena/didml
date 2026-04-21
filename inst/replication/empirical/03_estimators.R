# =============================================================================
# 03_estimators.R -- Stage 3: Compute all 5 estimators with 3 trimming variants
#
# Computes:
#   1. Wald         -- cell means, no covariates, no trimming
#   2. Wald_X       -- full-sample Lasso, no cross-fitting, no trimming
#   3. TC_X         -- full-sample Lasso, TC estimand, no cross-fitting, no trimming
#   4. DR-Wald      -- cross-fitted Lasso + FSIF, 3 trimming variants:
#                      (a) data-driven optimal (Remark 2, one-sided)
#                      (b) no-trim (alpha=0)
#                      (c) fixed symmetric [0.10, 0.90]
#   5. DR-TC        -- cross-fitted Lasso + FSIF, TC estimand, same 3 trims
#
# Total: 9 estimate columns (3 non-DR + 3 DR-Wald trims + 3 DR-TC trims)
#
# Also computes:
#   - State-clustered SEs for all estimators (Fix #3)
#   - Pre-trend / placebo tests (Fix #4):
#     (a) Placebo outcome (employment) Wald-DID
#     (b) Propensity-split falsification test
#
# Input:  data/clean/medicaid_analysis_v2.rds
#         data/clean/empirical_nuisance.rds
# Output: output/simulations/empirical_estimates.rds
#
# Author:  Andres Mena (Brown University)
# Project: fuzzyDDML
# =============================================================================

source("code/R/sim/00_config.R")
source("code/R/04_estimators.R")

# =============================================================================
# CLUSTERED VARIANCE (Fix #3: state-level treatment requires clustering)
# =============================================================================

# V_clustered = (1/N^2) * sum_s (sum_{i in s} psi_i)^2
# SE_clustered = sqrt(V_clustered / G_hat^2)
compute_clustered_se <- function(scores, cluster_ids, G_hat) {
  N <- length(scores)
  clusters <- unique(cluster_ids)
  cluster_sums <- tapply(scores, cluster_ids, sum)
  V_clust <- sum(cluster_sums^2) / N^2
  sqrt(V_clust / (G_hat^2))
}

# =============================================================================
# DATA-DRIVEN TRIMMING: Remark 2 (same function as in sim/03_estimators.R)
# =============================================================================
ALPHA_GRID_DD <- seq(0, 0.50, by = 0.01)

compute_data_driven_alpha <- function(pG_raw, DID_D_hat) {
  n <- length(pG_raw)

  # Estimate g_hat(e) = E[DID_D | e(X)=e] via local linear regression
  span_val <- max(0.3, min(0.75, 50 / n))
  g_fit <- tryCatch({
    loess(DID_D_hat ~ pG_raw, span = span_val, degree = 1)
  }, error = function(e) NULL)

  if (!is.null(g_fit)) {
    g_hat <- pmax(predict(g_fit, newdata = data.frame(pG_raw = pG_raw)), 0.001)
    g_hat[is.na(g_hat)] <- mean(DID_D_hat, na.rm = TRUE)
  } else {
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
# 1. LOAD DATA AND NUISANCE
# =============================================================================

message("=== 03_estimators.R: Compute 5 estimators x 3 trims ===\n")

args <- commandArgs(trailingOnly = TRUE)
use_duflo <- "--duflo" %in% args

if (use_duflo) {
  dat_file  <- file.path(DIR_EMPIRICAL, "duflo_analysis.rds")
  nuis_file <- file.path(DIR_EMPIRICAL, "empirical_nuisance_duflo.rds")
} else {
  dat_file  <- file.path(DIR_EMPIRICAL, "medicaid_analysis_v2.rds")
  nuis_file <- file.path(DIR_EMPIRICAL, "empirical_nuisance_medicaid.rds")
}
if (!file.exists(dat_file)) stop("Run 01_data.R first: ", dat_file)
if (!file.exists(nuis_file)) stop("Run 02_nuisance.R first: ", nuis_file)

dat  <- readRDS(dat_file)
nuis <- readRDS(nuis_file)

Y  <- dat$Y
D  <- dat$D
G  <- dat$G
Ti <- dat$Ti
N  <- length(Y)
pT <- nuis$pT
GT <- G * Ti
p_11 <- mean(GT)

# State identifiers for clustering (ACA is a state-level treatment)
STATEFIP <- if (!is.null(dat$STATEFIP)) dat$STATEFIP else dat$df$STATEFIP
n_states <- length(unique(STATEFIP))
message(sprintf("  Clustering: %d states", n_states))

fs <- nuis$fs
cf <- nuis$cf

message(sprintf("  N = %s, p = %d", format(N, big.mark = ","), dat$metadata$p))

# =============================================================================
# 2. WALD (cell means, no covariates)
# =============================================================================

message("\n1. Wald (cell means)...")

W_df <- data.frame(Y = Y, D = D, G = G, Ti = Ti)

# Cell means
mY_cell <- mD_cell <- list()
for (g_val in 0:1) for (t_val in 0:1) {
  k <- paste0(g_val, t_val)
  idx <- G == g_val & Ti == t_val
  mY_cell[[k]] <- rep(mean(Y[idx]), N)
  mD_cell[[k]] <- rep(mean(D[idx]), N)
}

DID_Y_cell <- Y - mY_cell[["10"]] - mY_cell[["01"]] + mY_cell[["00"]]
DID_D_cell <- D - mD_cell[["10"]] - mD_cell[["01"]] + mD_cell[["00"]]

wald_est <- sum(GT * DID_Y_cell) / sum(GT * DID_D_cell)

# SE via CLT on the cell-mean Wald ratio
# Use the efficient IF for the cell-mean case (degenerate propensity = Pr(G))
pG_cell <- rep(mean(G), N)
pi11_c <- pG_cell * pT; pi10_c <- pG_cell * (1 - pT)
pi01_c <- (1 - pG_cell) * pT; pi00_c <- (1 - pG_cell) * (1 - pT)
w10_c <- G * (1 - Ti) * pi11_c / (pi10_c * p_11)
w01_c <- (1 - G) * Ti * pi11_c / (pi01_c * p_11)
w00_c <- (1 - G) * (1 - Ti) * pi11_c / (pi00_c * p_11)

eps_c <- Y - wald_est * D
me10_c <- mY_cell[["10"]] - wald_est * mD_cell[["10"]]
me01_c <- mY_cell[["01"]] - wald_est * mD_cell[["01"]]
me00_c <- mY_cell[["00"]] - wald_est * mD_cell[["00"]]
psi_wald <- (GT / p_11) * (eps_c - me10_c - me01_c + me00_c) -
  w10_c * (eps_c - me10_c) - w01_c * (eps_c - me01_c) + w00_c * (eps_c - me00_c)
Gh_wald <- -mean((GT / p_11) * DID_D_cell -
  w10_c * (D - mD_cell[["10"]]) - w01_c * (D - mD_cell[["01"]]) + w00_c * (D - mD_cell[["00"]]))
se_wald <- sqrt(mean(psi_wald^2) / (Gh_wald^2) / N)

se_wald_clust <- compute_clustered_se(psi_wald, STATEFIP, Gh_wald)

result_wald <- list(estimate = wald_est, se = se_wald, se_clustered = se_wald_clust,
                    ci_lower = wald_est - 1.96 * se_wald,
                    ci_upper = wald_est + 1.96 * se_wald,
                    ci_lower_clust = wald_est - 1.96 * se_wald_clust,
                    ci_upper_clust = wald_est + 1.96 * se_wald_clust)

message(sprintf("  Wald: %.4f (SE %.4f, clustered SE %.4f)", wald_est, se_wald, se_wald_clust))

# =============================================================================
# 3. WALD_X (full-sample Lasso, no cross-fitting, no FSIF)
# =============================================================================

message("2. Wald_X (full-sample Lasso)...")

# Use full-sample nuisance for point estimate, efficient IF for SE
mY_fs <- list("10" = fs$m_Y_10, "01" = fs$m_Y_01, "00" = fs$m_Y_00)
mD_fs <- list("10" = fs$m_D_10, "01" = fs$m_D_01, "00" = fs$m_D_00)

DID_Y_fs <- Y - mY_fs[["10"]] - mY_fs[["01"]] + mY_fs[["00"]]
DID_D_fs <- D - mD_fs[["10"]] - mD_fs[["01"]] + mD_fs[["00"]]
waldx_est <- sum(GT * DID_Y_fs) / sum(GT * DID_D_fs)

# SE via efficient IF using full-sample pG (symmetric trimming at PG_TRIM_DEFAULT)
pG_fs_trim <- pmax(pmin(fs$pG_raw, 1 - PG_TRIM_DEFAULT), PG_TRIM_DEFAULT)
pi11_f <- pG_fs_trim * pT; pi10_f <- pG_fs_trim * (1 - pT)
pi01_f <- (1 - pG_fs_trim) * pT; pi00_f <- (1 - pG_fs_trim) * (1 - pT)
w10_f <- G * (1 - Ti) * pi11_f / (pi10_f * p_11)
w01_f <- (1 - G) * Ti * pi11_f / (pi01_f * p_11)
w00_f <- (1 - G) * (1 - Ti) * pi11_f / (pi00_f * p_11)

eps_f <- Y - waldx_est * D
me10_f <- mY_fs[["10"]] - waldx_est * mD_fs[["10"]]
me01_f <- mY_fs[["01"]] - waldx_est * mD_fs[["01"]]
me00_f <- mY_fs[["00"]] - waldx_est * mD_fs[["00"]]
psi_fx <- (GT / p_11) * (eps_f - me10_f - me01_f + me00_f) -
  w10_f * (eps_f - me10_f) - w01_f * (eps_f - me01_f) + w00_f * (eps_f - me00_f)
Gh_fx <- -mean((GT / p_11) * DID_D_fs -
  w10_f * (D - mD_fs[["10"]]) - w01_f * (D - mD_fs[["01"]]) + w00_f * (D - mD_fs[["00"]]))
se_waldx <- sqrt(mean(psi_fx^2) / (Gh_fx^2) / N)

se_waldx_clust <- compute_clustered_se(psi_fx, STATEFIP, Gh_fx)

result_waldx <- list(estimate = waldx_est, se = se_waldx, se_clustered = se_waldx_clust,
                     ci_lower = waldx_est - 1.96 * se_waldx,
                     ci_upper = waldx_est + 1.96 * se_waldx,
                     ci_lower_clust = waldx_est - 1.96 * se_waldx_clust,
                     ci_upper_clust = waldx_est + 1.96 * se_waldx_clust)

message(sprintf("  Wald_X: %.4f (SE %.4f, clustered SE %.4f)", waldx_est, se_waldx, se_waldx_clust))

# =============================================================================
# 4. TC_X (full-sample Lasso, TC estimand, no cross-fitting)
# =============================================================================

message("3. TC_X (full-sample Lasso, TC estimand)...")

delta_1_fs <- fs$mu_Y_101 - fs$mu_Y_100
delta_0_fs <- fs$mu_Y_001 - fs$mu_Y_000
gY_tc_fs <- Y - fs$m_Y_10 - delta_1_fs * fs$m_D_10 - (1 - fs$m_D_10) * delta_0_fs
gD_tc_fs <- D - fs$m_D_10

tcx_est <- sum(GT * gY_tc_fs) / sum(GT * gD_tc_fs)

# SE via the TC moment function
psi_tc_fs <- (GT / p_11) * (gY_tc_fs - gD_tc_fs * tcx_est)
Gh_tc_fs <- -mean((GT / p_11) * gD_tc_fs)
se_tcx <- sqrt(mean(psi_tc_fs^2) / (Gh_tc_fs^2) / N)

se_tcx_clust <- compute_clustered_se(psi_tc_fs, STATEFIP, Gh_tc_fs)

result_tcx <- list(estimate = tcx_est, se = se_tcx, se_clustered = se_tcx_clust,
                   ci_lower = tcx_est - 1.96 * se_tcx,
                   ci_upper = tcx_est + 1.96 * se_tcx,
                   ci_lower_clust = tcx_est - 1.96 * se_tcx_clust,
                   ci_upper_clust = tcx_est + 1.96 * se_tcx_clust)

message(sprintf("  TC_X: %.4f (SE %.4f, clustered SE %.4f)", tcx_est, se_tcx, se_tcx_clust))

# =============================================================================
# 5. DR-WALD (cross-fitted + FSIF, 3 trimming variants)
# =============================================================================

message("4. DR-Wald (3 trims)...")

# Compute data-driven alpha
DID_D_hat_cf <- D - cf$m_D_10 - cf$m_D_01 + cf$m_D_00
alpha_optimal <- compute_data_driven_alpha(cf$pG_raw, DID_D_hat_cf)
message(sprintf("  Data-driven alpha: %.3f", alpha_optimal))

# Compute the R_hat curve for later plotting
R_hat_curve <- data.frame(alpha = ALPHA_GRID_DD, R_hat = NA_real_)
g_fit_empirical <- tryCatch({
  loess(DID_D_hat_cf ~ cf$pG_raw, span = max(0.3, min(0.75, 50 / N)), degree = 1)
}, error = function(e) NULL)
g_hat_empirical <- if (!is.null(g_fit_empirical)) {
  g_tmp <- pmax(predict(g_fit_empirical, newdata = data.frame(pG_raw = cf$pG_raw)), 0.001)
  g_tmp[is.na(g_tmp)] <- mean(DID_D_hat_cf, na.rm = TRUE)
  g_tmp
} else {
  rep(mean(DID_D_hat_cf), N)
}

for (i in seq_along(ALPHA_GRID_DD)) {
  a <- ALPHA_GRID_DD[i]
  keep <- cf$pG_raw <= (1 - a)
  if (sum(keep) < 10) next
  num <- sum(cf$pG_raw[keep] / (1 - cf$pG_raw[keep]))
  den <- sum(cf$pG_raw[keep] * g_hat_empirical[keep])^2
  if (den > 1e-10) R_hat_curve$R_hat[i] <- num / den
}

# Three trimming approaches
trim_configs <- list(
  optimal = list(name = "optimal", alpha = alpha_optimal, one_sided = TRUE),
  notrim  = list(name = "no-trim", alpha = 0,             one_sided = TRUE),
  fixed   = list(name = "fixed",   alpha = PG_TRIM_DEFAULT, one_sided = FALSE)
)

result_drwald <- list()

for (tc_name in names(trim_configs)) {
  tc <- trim_configs[[tc_name]]

  # Apply trimming to pG
  if (tc$one_sided) {
    pG_trimmed <- pmin(cf$pG_raw, 1 - tc$alpha)
    pG_trimmed <- pmax(pG_trimmed, 0.001)
  } else {
    pG_trimmed <- pmax(pmin(cf$pG_raw, 1 - tc$alpha), tc$alpha)
  }

  # Reconstruct cell probabilities
  pred_w <- cf
  pred_w$pi_11 <- pG_trimmed * pT
  pred_w$pi_10 <- pG_trimmed * (1 - pT)
  pred_w$pi_01 <- (1 - pG_trimmed) * pT
  pred_w$pi_00 <- (1 - pG_trimmed) * (1 - pT)

  res <- tryCatch(compute_dd_wald(W_df, pred_w), error = function(e) {
    list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
         scores = rep(NA, N), variance = NA)
  })

  # Clustered SE using individual scores from compute_dd_wald
  if (!is.null(res$scores) && !all(is.na(res$scores))) {
    G_hat_w <- -mean((GT / p_11) * DID_D_hat_cf -
      G * (1 - Ti) * pred_w$pi_11 / (pred_w$pi_10 * p_11) * (D - cf$m_D_10) -
      (1 - G) * Ti * pred_w$pi_11 / (pred_w$pi_01 * p_11) * (D - cf$m_D_01) +
      (1 - G) * (1 - Ti) * pred_w$pi_11 / (pred_w$pi_00 * p_11) * (D - cf$m_D_00))
    res$se_clustered <- compute_clustered_se(res$scores, STATEFIP, G_hat_w)
    res$ci_lower_clust <- res$estimate - 1.96 * res$se_clustered
    res$ci_upper_clust <- res$estimate + 1.96 * res$se_clustered
  } else {
    res$se_clustered <- NA
    res$ci_lower_clust <- NA
    res$ci_upper_clust <- NA
  }

  result_drwald[[tc_name]] <- res
  message(sprintf("  DR-Wald (%s, alpha=%.2f): %.4f (SE %.4f, clustered SE %.4f)",
                  tc$name, tc$alpha, res$estimate, res$se,
                  ifelse(is.na(res$se_clustered), NA, res$se_clustered)))
}

# =============================================================================
# 6. DR-TC (cross-fitted + FSIF, TC estimand, 3 trimming variants)
# =============================================================================

message("5. DR-TC (3 trims)...")

result_drtc <- list()

for (tc_name in names(trim_configs)) {
  tc <- trim_configs[[tc_name]]

  if (tc$one_sided) {
    pG_trimmed <- pmin(cf$pG_raw, 1 - tc$alpha)
    pG_trimmed <- pmax(pG_trimmed, 0.001)
  } else {
    pG_trimmed <- pmax(pmin(cf$pG_raw, 1 - tc$alpha), tc$alpha)
  }

  pred_t <- cf
  pred_t$pi_11 <- pG_trimmed * pT
  pred_t$pi_10 <- pG_trimmed * (1 - pT)
  pred_t$pi_01 <- (1 - pG_trimmed) * pT
  pred_t$pi_00 <- (1 - pG_trimmed) * (1 - pT)

  res <- tryCatch(compute_dd_tc(W_df, pred_t), error = function(e) {
    list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
         scores = rep(NA, N), variance = NA)
  })

  # Clustered SE using individual scores from compute_dd_tc
  if (!is.null(res$scores) && !all(is.na(res$scores))) {
    # G_hat for TC: -mean(denominator_i) from compute_dd_tc
    DID_D_tc <- D - cf$m_D_10
    G_hat_t <- -mean((GT / p_11) * DID_D_tc -
      G * (1 - Ti) * pred_t$pi_11 / (pred_t$pi_10 * p_11) * (D - cf$m_D_10))
    res$se_clustered <- compute_clustered_se(res$scores, STATEFIP, G_hat_t)
    res$ci_lower_clust <- res$estimate - 1.96 * res$se_clustered
    res$ci_upper_clust <- res$estimate + 1.96 * res$se_clustered
  } else {
    res$se_clustered <- NA
    res$ci_lower_clust <- NA
    res$ci_upper_clust <- NA
  }

  result_drtc[[tc_name]] <- res
  message(sprintf("  DR-TC (%s, alpha=%.2f): %.4f (SE %.4f, clustered SE %.4f)",
                  tc$name, tc$alpha, res$estimate, res$se,
                  ifelse(is.na(res$se_clustered), NA, res$se_clustered)))
}

# =============================================================================
# 7. COMPUTE SUMMARY STATISTICS FOR TABLE
# =============================================================================

did_Y <- dat$metadata$did_Y
did_D <- dat$metadata$did_D

# =============================================================================
# 7b. PRE-TREND / PLACEBO TEST (Fix #4)
# =============================================================================
#
# With only 2013 and 2015 available, we cannot compute a placebo DID using
# two pre-periods. Instead, we run two falsification tests:
#
# (a) Placebo outcome: estimate the Wald-DID on employment (not targeted by
#     ACA Medicaid expansion). Under no spillovers, this should be zero.
#
# (b) Control-group falsification: within G=0 states, split at the median
#     propensity score and compute a "pseudo-DID" treating high-propensity
#     control units as a placebo treatment group. Under parallel trends,
#     the reduced-form DID_Y should be zero.

message("\n7b. Pre-trend / placebo tests...")

# --- (a) Placebo outcome: employment ---
placebo_Y <- as.numeric(dat$df$EMPSTAT == 1)  # employed indicator
DID_empY <- mean(placebo_Y[G == 1 & Ti == 1]) - mean(placebo_Y[G == 1 & Ti == 0]) -
           (mean(placebo_Y[G == 0 & Ti == 1]) - mean(placebo_Y[G == 0 & Ti == 0]))
DID_empD <- did_D  # same first stage

placebo_wald <- DID_empY / DID_empD

# SE for placebo Wald via delta method on cell means
n_cells <- table(G, Ti)
var_empY <- var(placebo_Y)
# Approximate SE using cell-mean variance
se_num <- sqrt(var_empY * sum(1 / n_cells))
se_placebo <- abs(se_num / DID_empD)

result_placebo_employment <- list(
  estimate = placebo_wald,
  se = se_placebo,
  did_empY = DID_empY,
  pvalue = 2 * pnorm(-abs(placebo_wald / se_placebo)),
  description = "Placebo outcome: employment (should be ~0 under no spillovers)"
)

message(sprintf("  Placebo (employment): %.4f (SE %.4f, p=%.3f)",
                placebo_wald, se_placebo, result_placebo_employment$pvalue))

# --- (b) Control-group falsification ---
# Among G=0 states, define pseudo-treatment by median propensity score
pG_median <- median(cf$pG_raw)
pseudo_G <- as.numeric(cf$pG_raw > pG_median)

# Compute pseudo-DID_Y and pseudo-DID_D using the pseudo-group indicator
pseudo_DID_Y <- mean(Y[pseudo_G == 1 & Ti == 1]) - mean(Y[pseudo_G == 1 & Ti == 0]) -
               (mean(Y[pseudo_G == 0 & Ti == 1]) - mean(Y[pseudo_G == 0 & Ti == 0]))

# Under parallel trends, the pseudo reduced-form DID should be zero
# Test H0: pseudo_DID_Y = 0
n_pseudo_cells <- c(sum(pseudo_G == 1 & Ti == 1), sum(pseudo_G == 1 & Ti == 0),
                    sum(pseudo_G == 0 & Ti == 1), sum(pseudo_G == 0 & Ti == 0))
se_pseudo <- sqrt(var(Y) * sum(1 / n_pseudo_cells))

result_placebo_propensity <- list(
  estimate = pseudo_DID_Y,
  se = se_pseudo,
  pvalue = 2 * pnorm(-abs(pseudo_DID_Y / se_pseudo)),
  pG_median = pG_median,
  description = "Propensity-split falsification: pseudo-DID_Y among all obs split at median e(X)"
)

message(sprintf("  Placebo (propensity split): DID_Y=%.4f (SE %.4f, p=%.3f)",
                pseudo_DID_Y, se_pseudo, result_placebo_propensity$pvalue))

# =============================================================================
# 8. SAVE
# =============================================================================

out <- list(
  wald      = result_wald,
  waldx     = result_waldx,
  tcx       = result_tcx,
  drwald    = result_drwald,
  drtc      = result_drtc,
  alpha_optimal = alpha_optimal,
  R_hat_curve = R_hat_curve,
  g_hat     = g_hat_empirical,
  pG_raw_cf = cf$pG_raw,
  did_Y     = did_Y,
  did_D     = did_D,
  N         = N,
  p         = dat$metadata$p,
  n_states  = n_states,
  trim_configs = trim_configs,
  placebo_employment = result_placebo_employment,
  placebo_propensity = result_placebo_propensity,
  timestamp = Sys.time()
)

outfile <- sprintf("output/simulations/empirical_estimates_%s.rds",
                   if (use_duflo) "duflo" else "medicaid")
saveRDS(out, outfile)
message(sprintf("\nSaved: %s", outfile))

# =============================================================================
# 9. PRINT SUMMARY
# =============================================================================

message("\n========== RESULTS SUMMARY ==========")
message(sprintf("  Wald:            %.4f  SE=%.4f  ClustSE=%.4f",
                result_wald$estimate, result_wald$se, result_wald$se_clustered))
message(sprintf("  Wald_X:          %.4f  SE=%.4f  ClustSE=%.4f",
                result_waldx$estimate, result_waldx$se, result_waldx$se_clustered))
message(sprintf("  TC_X:            %.4f  SE=%.4f  ClustSE=%.4f",
                result_tcx$estimate, result_tcx$se, result_tcx$se_clustered))
for (tc_name in names(result_drwald)) {
  r <- result_drwald[[tc_name]]
  message(sprintf("  DR-Wald(%s): %.4f  SE=%.4f  ClustSE=%.4f",
                  tc_name, r$estimate, r$se,
                  ifelse(is.na(r$se_clustered), NA, r$se_clustered)))
}
for (tc_name in names(result_drtc)) {
  r <- result_drtc[[tc_name]]
  message(sprintf("  DR-TC(%s):   %.4f  SE=%.4f  ClustSE=%.4f",
                  tc_name, r$estimate, r$se,
                  ifelse(is.na(r$se_clustered), NA, r$se_clustered)))
}
message(sprintf("\n  DID_Y = %.4f, DID_D = %.4f", did_Y, did_D))
message(sprintf("  N = %s, p = %d, States = %d", format(N, big.mark = ","), dat$metadata$p, n_states))
message(sprintf("  Optimal alpha = %.3f", alpha_optimal))

message("\n=== Stage 3 complete ===")
