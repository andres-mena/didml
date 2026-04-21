# =============================================================================
# 01_dgp.R — Stage 1: Generate DGP data from pre-registered seeds
#
# Generates M replications of DGP v6a at each sample size, augments X to p=100.
# NO estimation. Data only.
#
# Input:  00_config.R (seeds, DGP params)
# Output: data/simulations/dgp/dgp_n{N}.rds
#
#   Structure per file:
#   $data     — list of M data.frames, each n x (9 + p) columns:
#               Y, Y0, Y1, D, Ti, G, V, g0, h, X1, ..., X{p}
#   $seeds    — integer vector of M pre-registered seeds
#   $params   — list of DGP parameters (tau, p_orig, p_aug, rho, alpha_corr)
#   $n, $M, $timestamp
#
# Usage:
#   Rscript code/R/sim/01_dgp.R              # all N_VALUES
#   Rscript code/R/sim/01_dgp.R --n 500      # single n
#
# =============================================================================

source("code/R/sim/00_config.R")
library(MASS)

args <- commandArgs(trailingOnly = TRUE)
n_filter <- NULL
for (i in seq_along(args)) {
  if (args[i] == "--n" && i < length(args)) n_filter <- as.integer(args[i + 1])
}
n_run <- if (!is.null(n_filter)) n_filter else N_VALUES

# =============================================================================
# DGP v6a — nonlinear DID with fuzzy compliance
# =============================================================================
gen_dgp_v6a <- function(n, tau = TAU, p = P_ORIG, rho = RHO) {
  Sigma <- toeplitz(rho^(0:(p - 1))) * 2.25
  X <- mvrnorm(n, rep(0, p), Sigma)
  colnames(X) <- paste0("X", 1:p)

  Ti <- rbinom(n, 1, 0.5)
  p_g <- plogis(0.6 * X[, 1] + 0.3 * X[, 3])
  G <- rbinom(n, 1, p_g)

  h_x <- 0.4 * X[, 1] + 0.2 * X[, 1]^2 + 0.15 * X[, 3]
  V <- G + 0.4 * X[, 1] + 0.2 * X[, 3] + rnorm(n)

  D0 <- Ti * (V > 2) + (1 - Ti) * (V > 2)
  D1 <- Ti * (V > 0) + (1 - Ti) * (V > 2)
  D <- G * D1 + (1 - G) * D0

  g0 <- plogis(X[, 1]) + 0.3 * X[, 3] + 0.3 * X[, 1] * X[, 3] +
    0.2 * X[, 2]^2 + 0.15 * X[, 4]

  Y0 <- 0.3 * G + h_x * Ti + g0 + rnorm(n)
  Y1 <- Y0 + tau + rnorm(n, 0, 0.5)
  Y <- D * Y1 + (1 - D) * Y0

  data.frame(Y = Y, Y0 = Y0, Y1 = Y1, D = D, Ti = Ti, G = G,
             V = V, g0 = g0, h = h_x, X)
}

# =============================================================================
# Augment X: add X21-X100 correlated with signal covariates X1-X5
# =============================================================================
augment_covariates <- function(W, p_target = P_AUG, alpha = ALPHA_CORR, aug_seed) {
  set.seed(aug_seed)
  n <- nrow(W)
  p_current <- sum(grepl("^X", names(W)))
  p_new <- p_target - p_current
  if (p_new <= 0) return(W)

  parents <- ((1:p_new - 1) %% 5) + 1
  sd_x <- sqrt(2.25)

  X_new <- matrix(NA_real_, n, p_new)
  for (j in 1:p_new) {
    eps <- rnorm(n, 0, sd_x)
    X_new[, j] <- alpha * W[[paste0("X", parents[j])]] + sqrt(1 - alpha^2) * eps
  }
  colnames(X_new) <- paste0("X", (p_current + 1):p_target)
  cbind(W, X_new)
}

# =============================================================================
# MAIN
# =============================================================================
cat("=============================================================================\n")
cat("01_dgp.R — Generate DGP data\n")
cat(sprintf("  n_values: %s | M=%d | p=%d→%d\n",
            paste(n_run, collapse = ","), M, P_ORIG, P_AUG))
cat("=============================================================================\n\n")

for (n_val in n_run) {
  seeds <- get_seeds(n_val)
  stopifnot(length(seeds) >= M)
  seeds <- seeds[1:M]

  cat(sprintf("[%s] n=%d: generating %d reps...\n",
              format(Sys.time(), "%H:%M:%S"), n_val, M))
  t0 <- Sys.time()

  data_list <- vector("list", M)
  for (r in 1:M) {
    set.seed(seeds[r])
    W <- gen_dgp_v6a(n = n_val)
    W <- augment_covariates(W, p_target = P_AUG, alpha = ALPHA_CORR,
                            aug_seed = AUG_OFFSET + seeds[r])
    data_list[[r]] <- W
  }

  outfile <- dgp_file(n_val)
  saveRDS(list(
    data   = data_list,
    seeds  = seeds,
    params = list(tau = TAU, p_orig = P_ORIG, p_aug = P_AUG,
                  rho = RHO, alpha_corr = ALPHA_CORR, dgp = "v6a"),
    n = n_val, M = M,
    timestamp = Sys.time()
  ), outfile)

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  mb <- round(file.size(outfile) / 1e6, 1)
  cat(sprintf("  Done in %.1fs → %s (%.1f MB)\n\n", elapsed, outfile, mb))
}

cat("=== Stage 1 complete ===\n")
