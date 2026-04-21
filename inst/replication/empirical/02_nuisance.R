# =============================================================================
# 02_nuisance.R -- Stage 2: Cross-fit all nuisance functions (EXPENSIVE)
#
# Estimates all nuisance components for the 5 empirical estimators:
#   - Conditional expectations m^R_{gt}(X) for R in {Y, D}, (g,t) != (1,1)
#   - Propensity scores Pr(G=1|X) via logit-Lasso
#   - TC treatment-conditional expectations mu^Y_{dgt}(X) in control group
#   - TC cell propensities pi_{dgt}(X)
#
# Produces BOTH full-sample (for Wald_X, TC_X) and cross-fitted (for DR-Wald,
# DR-TC) predictions. Saves raw pG (untrimmed) so trimming can be applied
# downstream in Stage 3.
#
# Input:  data/clean/medicaid_analysis_v2.rds
# Output: data/clean/empirical_nuisance.rds
#
# Author:  Andres Mena (Brown University)
# Project: fuzzyDDML
# =============================================================================

source("code/R/sim/00_config.R")

library(glmnet)

set.seed(EMPIRICAL_SEED)

K <- EMPIRICAL_K

# =============================================================================
# 1. LOAD DATA
# =============================================================================

message("=== 02_nuisance.R: Cross-fit nuisance functions ===\n")

# Select dataset: --duflo flag uses SUPAS data, default uses Medicaid/ACS
args <- commandArgs(trailingOnly = TRUE)
use_duflo <- "--duflo" %in% args

if (use_duflo) {
  infile <- file.path(DIR_EMPIRICAL, "duflo_analysis.rds")
  outfile_suffix <- "duflo"
} else {
  infile <- file.path(DIR_EMPIRICAL, "medicaid_analysis_v2.rds")
  outfile_suffix <- "medicaid"
}
if (!file.exists(infile)) stop("Run 01_data first: ", infile)

dat <- readRDS(infile)
Y  <- dat$Y
D  <- dat$D
G  <- dat$G
Ti <- dat$Ti
X  <- dat$X_matrix
STATEFIP <- dat$STATEFIP
N  <- length(Y)
p  <- ncol(X)
pT <- mean(Ti)

# Cluster-level folds are theoretically preferred under cluster dependence
# (Chiang et al. 2022, JBES §3.1.2): they preserve cross-fit independence
# between \hat{\gamma}_l and I_l, which drives the C^* term to zero.
#
# In the INPRES application (253 districts, p=197 district-invariant
# covariates) cluster folds nevertheless collapse the Lasso nuisance in
# finite samples. Empirical diagnostic (see explorations/cluster_folds_
# comparison/, 2026-04-20): outcome pseudo-R^2 = 0.000 for m_Y_{10}, m_Y_{01},
# m_Y_{00}; cross-fitted propensities on [0, 1] boundary; DR-Wald blows up
# to ~90 from ~0.14. Mechanism: cluster-blocked inner CV sees too few
# distinct X-vectors per training fold relative to p, so lambda -> infty
# and the Lasso shrinks to the intercept. Observation-level folds recover
# informative nuisance fits (pseudo-R^2 > 0, bounded propensities).
#
# Default therefore uses observation-level folds; --cluster-folds is
# opt-in for settings where p << C or X has within-cluster variation.
use_cluster_folds <- ("--cluster-folds" %in% args) && !is.null(STATEFIP)
if (use_cluster_folds) {
  C <- length(unique(STATEFIP))
  message(sprintf("  N = %s, p = %d, K = %d, C = %d (cluster-level folds)",
                  format(N, big.mark = ","), p, K, C))
} else {
  message(sprintf("  N = %s, p = %d, K = %d (observation-level folds)",
                  format(N, big.mark = ","), p, K))
}

# =============================================================================
# 2. FULL-SAMPLE NUISANCE (for Wald_X, TC_X -- no cross-fitting)
# =============================================================================

message("\nPhase 1: Full-sample Lasso (Wald_X, TC_X)...")
t0 <- Sys.time()

fs <- list()

# --- m^R_{gt}(X) = E[R | G=g, T=t, X] via Lasso ---
for (R_name in c("Y", "D")) {
  R_vec <- if (R_name == "Y") Y else D
  for (g_val in 0:1) {
    for (t_val in 0:1) {
      if (g_val == 1 && t_val == 1) next
      key <- sprintf("m_%s_%d%d", R_name, g_val, t_val)
      idx <- which(G == g_val & Ti == t_val)

      if (length(idx) < 20) {
        fs[[key]] <- rep(mean(R_vec[idx]), N)
        next
      }

      fs[[key]] <- tryCatch({
        mod <- cv.glmnet(X[idx, , drop = FALSE], R_vec[idx], alpha = 1, nfolds = 5)
        as.vector(predict(mod, newx = X, s = "lambda.min"))
      }, error = function(e) rep(mean(R_vec[idx]), N))
    }
  }
}

# --- Full-sample Pr(G=1|X) via logit-Lasso ---
fs$pG_raw <- tryCatch({
  mod <- cv.glmnet(X, G, alpha = 1, family = "binomial", nfolds = 5)
  as.vector(predict(mod, newx = X, s = "lambda.min", type = "response"))
}, error = function(e) rep(mean(G), N))

# --- TC nuisance: mu^Y_{dgt}(X) = E[Y | D=d, G=0, T=t, X] ---
for (d_val in 0:1) {
  for (t_val in 0:1) {
    key <- sprintf("mu_Y_%d0%d", d_val, t_val)
    idx <- which(D == d_val & G == 0 & Ti == t_val)

    if (length(idx) < 20) {
      fs[[key]] <- rep(if (length(idx) > 0) mean(Y[idx]) else 0, N)
      next
    }

    fs[[key]] <- tryCatch({
      mod <- cv.glmnet(X[idx, , drop = FALSE], Y[idx], alpha = 1, nfolds = 5)
      as.vector(predict(mod, newx = X, s = "lambda.min"))
    }, error = function(e) rep(mean(Y[idx]), N))
  }
}

t_fs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
message(sprintf("  Full-sample done in %.1fs", t_fs))

# =============================================================================
# 3. CROSS-FITTED NUISANCE (for DR-Wald, DR-TC)
# =============================================================================

message(sprintf("\nPhase 2: %d-fold cross-fitting (DR-Wald, DR-TC)...", K))
t1 <- Sys.time()

if (use_cluster_folds) {
  cluster_ids   <- unique(STATEFIP)
  cluster_fold  <- sample(rep(1:K, length.out = length(cluster_ids)))
  fold_ids      <- cluster_fold[match(STATEFIP, cluster_ids)]
  message(sprintf("  Cluster-level folds: %d clusters split into %d folds",
                  length(cluster_ids), K))
} else {
  fold_ids <- sample(rep(1:K, length.out = N))
}

# Initialize storage
cf <- list(fold_ids = fold_ids)
cf_fields <- c("m_Y_10", "m_Y_01", "m_Y_00",
               "m_D_10", "m_D_01", "m_D_00",
               "pG_raw",
               "mu_Y_101", "mu_Y_100", "mu_Y_001", "mu_Y_000",
               "pi_101", "pi_100", "pi_001", "pi_000")
for (nm in cf_fields) cf[[nm]] <- numeric(N)

for (k in 1:K) {
  tr <- which(fold_ids != k)
  va <- which(fold_ids == k)
  X_tr <- X[tr, , drop = FALSE]
  X_va <- X[va, , drop = FALSE]
  n_va <- length(va)

  message(sprintf("  Fold %d/%d (train=%d, valid=%d)...", k, K, length(tr), n_va))

  # --- Conditional expectations m^R_{gt}(X) ---
  for (R_name in c("Y", "D")) {
    R_vec <- if (R_name == "Y") Y else D
    for (g_val in 0:1) {
      for (t_val in 0:1) {
        if (g_val == 1 && t_val == 1) next
        key <- sprintf("m_%s_%d%d", R_name, g_val, t_val)
        si <- which(G[tr] == g_val & Ti[tr] == t_val)

        if (length(si) < 20) {
          cf[[key]][va] <- if (length(si) > 0) mean(R_vec[tr[si]]) else 0
          next
        }

        cf[[key]][va] <- tryCatch({
          mod <- cv.glmnet(X_tr[si, , drop = FALSE], R_vec[tr[si]],
                           alpha = 1, nfolds = 5)
          as.vector(predict(mod, newx = X_va, s = "lambda.min"))
        }, error = function(e) rep(mean(R_vec[tr[si]]), n_va))
      }
    }
  }

  # --- Propensity score Pr(G=1|X) ---
  cf$pG_raw[va] <- tryCatch({
    mod <- cv.glmnet(X_tr, G[tr], alpha = 1, family = "binomial", nfolds = 5)
    as.vector(predict(mod, newx = X_va, s = "lambda.min", type = "response"))
  }, error = function(e) rep(mean(G[tr]), n_va))

  # --- TC nuisance: mu^Y_{dgt}(X) ---
  for (d_val in 0:1) {
    for (t_val in 0:1) {
      key_mu <- sprintf("mu_Y_%d0%d", d_val, t_val)
      si <- which(D[tr] == d_val & G[tr] == 0 & Ti[tr] == t_val)

      if (length(si) < 20) {
        cf[[key_mu]][va] <- if (length(si) > 0) mean(Y[tr[si]]) else 0
        next
      }

      cf[[key_mu]][va] <- tryCatch({
        mod <- cv.glmnet(X_tr[si, , drop = FALSE], Y[tr[si]],
                         alpha = 1, nfolds = 5)
        as.vector(predict(mod, newx = X_va, s = "lambda.min"))
      }, error = function(e) rep(mean(Y[tr[si]]), n_va))
    }
  }

  # --- TC cell propensities: pi_{dgt}(X) = Pr(D=d, G=0, T=t | X) ---
  for (d_val in 0:1) {
    for (t_val in 0:1) {
      key_pi <- sprintf("pi_%d0%d", d_val, t_val)
      L_ind <- as.numeric(D[tr] == d_val & G[tr] == 0 & Ti[tr] == t_val)

      cf[[key_pi]][va] <- tryCatch({
        mod <- cv.glmnet(X_tr, L_ind, alpha = 1, family = "binomial", nfolds = 5)
        p_hat <- as.vector(predict(mod, newx = X_va, s = "lambda.min", type = "response"))
        pmax(pmin(p_hat, 0.99), 0.01)
      }, error = function(e) rep(pmax(mean(L_ind), 0.01), n_va))
    }
  }
}

# Derive cell probabilities from pG (untrimmed for now -- trimming in Stage 3)
cf$pi_11 <- cf$pG_raw * pT
cf$pi_10 <- cf$pG_raw * (1 - pT)
cf$pi_01 <- (1 - cf$pG_raw) * pT
cf$pi_00 <- (1 - cf$pG_raw) * (1 - pT)

t_cf <- round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1)
message(sprintf("  Cross-fitting done in %.1f min", t_cf))

# =============================================================================
# 4. DIAGNOSTICS
# =============================================================================

message("\n--- Nuisance diagnostics ---")
message(sprintf("  pG_raw (cross-fitted): mean=%.3f, sd=%.3f, range=[%.3f, %.3f]",
                mean(cf$pG_raw), sd(cf$pG_raw), min(cf$pG_raw), max(cf$pG_raw)))
message(sprintf("  pG_raw (full-sample):  mean=%.3f, sd=%.3f, range=[%.3f, %.3f]",
                mean(fs$pG_raw), sd(fs$pG_raw), min(fs$pG_raw), max(fs$pG_raw)))

# Check m^Y_{gt} quality: in-sample R^2 proxy
for (gt in c("10", "01", "00")) {
  key <- paste0("m_Y_", gt)
  g_val <- as.integer(substr(gt, 1, 1))
  t_val <- as.integer(substr(gt, 2, 2))
  idx <- G == g_val & Ti == t_val
  r2 <- 1 - mean((Y[idx] - cf[[key]][idx])^2) / var(Y[idx])
  message(sprintf("  m_Y_%s cross-fitted pseudo-R2: %.3f", gt, max(r2, 0)))
}

# =============================================================================
# 5. SAVE
# =============================================================================

out <- list(
  fs        = fs,
  cf        = cf,
  N         = N,
  p         = p,
  K         = K,
  pT        = pT,
  timestamp = Sys.time()
)

outfile <- file.path(DIR_EMPIRICAL, sprintf("empirical_nuisance_%s.rds", outfile_suffix))
saveRDS(out, outfile)
message(sprintf("\nSaved: %s (%.1f MB)", outfile, file.size(outfile) / 1e6))

message("\n=== Stage 2 complete ===")
