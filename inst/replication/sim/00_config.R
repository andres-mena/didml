# =============================================================================
# 00_config.R — Single source of truth for the simulation pipeline
#
# All parameters, seeds, paths, and estimator definitions live here.
# Every other script in code/R/sim/ sources this file.
#
# Author:  Andres Mena (Brown University)
# Project: fuzzyDDML
# =============================================================================

# --- Sample sizes and replications ---
N_VALUES  <- c(100, 500, 1000, 2000, 5000)
M         <- 1000L
TAU       <- 1

# --- DGP parameters ---
P_ORIG    <- 20L      # original covariates
P_AUG     <- 100L     # augmented covariates (X21-X100 correlated with X1-X5)
ALPHA_CORR <- 0.3     # correlation of augmented with parent signal covariates
RHO       <- 0.5      # Toeplitz correlation for X1-X20

# --- ML hyperparameters ---
K_FOLDS   <- 5L       # cross-fitting folds
P_SIEVE   <- 20L      # sieve basis uses X1-X20 (degree-2 polynomial = 231 terms)

# --- Alternative ML methods (pre-registered, 2026-03-25) ---
# Random Forest
RF_NTREE  <- 500L
RF_MTRY   <- NULL      # defaults to floor(p/3) for regression, floor(sqrt(p)) for classification
RF_NODESIZE <- 5L

# Neural Network (nnet: single hidden layer)
NN_SIZE   <- 5L        # hidden units
NN_DECAY  <- 0.01      # weight decay (L2 regularization)
NN_MAXIT  <- 200L
NN_LINOUT <- TRUE      # linear output for regression; FALSE for classification

# --- Trimming ---
PG_TRIM_DEFAULT <- 0.10  # default propensity score clipping
TC_TRIM         <- 0.01  # TC propensity floor
TRIM_GRID       <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25)

# --- Data-driven trimming (Corollary 2.1, Remark 2) ---
# DR-Wald and DR-TC use per-replication data-driven alpha computed in 03_estimators.R.
# The formula is R_hat(alpha) with estimated g_hat(e) = E[DID_D | e(X)=e] via loess.
# One-sided: pG <- min(pG, 1 - alpha). No lower bound needed (Corollary 2.1).
# Mean alpha by n (DGP v6a, M=1000, 2026-03-25):
#   n=100: 0.050 | n=500: 0.028 | n=1000: 0.038 | n=2000: 0.071 | n=5000: 0.073
# Non-DR estimators use PG_TRIM_DEFAULT (symmetric).

# --- Computation ---
N_CORES <- {
  env_val <- Sys.getenv("SIM_CORES", unset = "")
  if (nzchar(env_val)) as.integer(env_val) else max(1, parallel::detectCores() - 2L)
}

# --- Pre-registered seeds ---
# These must NEVER change. They reproduce the exact DGP draws.
ORIG_BASE     <- 20260324L
ORIG_M        <- 500L
ORIG_N_VALUES <- c(50L, 100L, 500L, 1000L)

EXT_A_BASE    <- 20270000L
EXT_A_STRIDE  <- 1000L
EXT_A_M       <- 500L

EXT_B_BASE    <- 20280000L
EXT_B_STRIDE  <- 1000L
EXT_B_M       <- 1000L
EXT_B_N_VALUES <- c(2000L, 5000L, 10000L)

AUG_OFFSET    <- 99000000L  # non-overlapping with all DGP seeds

get_seeds <- function(n) {
  if (n %in% ORIG_N_VALUES) {
    n_idx <- match(n, ORIG_N_VALUES)
    seeds_orig <- ORIG_BASE + (n_idx - 1L) * ORIG_M + seq_len(ORIG_M)
    seeds_extA <- EXT_A_BASE + (n_idx - 1L) * EXT_A_STRIDE + seq_len(EXT_A_M)
    return(c(seeds_orig, seeds_extA))
  } else if (n %in% EXT_B_N_VALUES) {
    n_idx <- match(n, EXT_B_N_VALUES)
    return(EXT_B_BASE + (n_idx - 1L) * EXT_B_STRIDE + seq_len(EXT_B_M))
  } else {
    stop(sprintf("Unknown n=%d — add to seed registry in 00_config.R", n))
  }
}

# --- Estimator registry ---
# Canonical names used throughout the pipeline
ESTIMATOR_NAMES <- c("wald", "waldx_lasso", "waldx_sieve",
                     "tcx_lasso", "tcx_sieve", "drwald", "drtc")

ESTIMATOR_LABELS <- c("Wald", "Wald$_X$(Lasso)", "Wald$_X$(Sieve)",
                      "TC$_X$(Lasso)", "TC$_X$(Sieve)", "DML-Wald", "DML-TC")

# Table 1 displays 5 estimators in this order (Sieve variants excluded)
TABLE1_NAMES  <- c("wald", "waldx_lasso", "drwald", "tcx_lasso", "drtc")
TABLE1_LABELS <- c("Wald", "Wald$_X$", "DML-Wald", "TC$_X$", "DML-TC")

# --- Directory structure ---
DIR_DGP      <- "data/simulations/dgp"
DIR_NUISANCE <- "data/simulations/nuisance"
DIR_ESTIMATES <- "output/simulations/estimates"
DIR_TABLES   <- "output/tables"
DIR_FIGURES  <- "output/figures"

for (d in c(DIR_DGP, DIR_NUISANCE, DIR_ESTIMATES, DIR_TABLES, DIR_FIGURES)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Empirical config ---
EMPIRICAL_SEED  <- 20260325L
EMPIRICAL_N_SUB <- 50000L
EMPIRICAL_K     <- 5L
DIR_EMPIRICAL   <- "data/clean"
dir.create(DIR_EMPIRICAL, recursive = TRUE, showWarnings = FALSE)

# --- File naming conventions ---
dgp_file       <- function(n) file.path(DIR_DGP, sprintf("dgp_n%d.rds", n))
nuisance_file  <- function(n, method = "lasso") {
  if (method == "lasso") file.path(DIR_NUISANCE, sprintf("nuisance_n%d.rds", n))
  else file.path(DIR_NUISANCE, sprintf("nuisance_%s_n%d.rds", method, n))
}
estimates_file <- function(n, method = "lasso") {
  if (method == "lasso") file.path(DIR_ESTIMATES, sprintf("estimates_n%d.rds", n))
  else file.path(DIR_ESTIMATES, sprintf("estimates_%s_n%d.rds", method, n))
}

# --- Sieve basis builder ---
build_sieve <- function(X_mat, p_s = P_SIEVE) {
  X_sub <- X_mat[, 1:min(p_s, ncol(X_mat)), drop = FALSE]
  p <- ncol(X_sub)
  basis <- X_sub
  for (j in 1:p) for (k in j:p) basis <- cbind(basis, X_sub[, j] * X_sub[, k])
  basis
}

cat(sprintf("00_config.R loaded: N=%s, M=%d, p=%d, K=%d, cores=%d\n",
            paste(N_VALUES, collapse = ","), M, P_AUG, K_FOLDS, N_CORES))
