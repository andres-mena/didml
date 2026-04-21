# =============================================================================
# 02_nuisance.R — Stage 2: Estimate nuisance functions (EXPENSIVE, run once)
#
# For each replication: fits Lasso, Sieve, and cross-fitted ML predictions.
# Saves raw data vectors (Y, D, G, Ti, V, h, g0, X) alongside predictions
# so that Stage 3 never needs to touch the multi-GB DGP files.
#
# Input:  data/simulations/dgp/dgp_n{N}.rds
# Output: data/simulations/nuisance/nuisance_n{N}.rds
#
#   Structure per file:
#   $nuisance — list of M lists, each containing:
#       ## Raw data (for estimator computation in Stage 3)
#       Y, D, G, Ti           — outcome/treatment/group/time vectors (length n)
#       V, h, g0              — latent variables for diagnostics
#       X                     — covariate matrix (n x p)
#       ## Full-sample predictions (for Wald_X, TC_X estimators)
#       fs_lasso_{Y,D}_{gt}   — Lasso predictions in each (g,t) cell
#       fs_sieve_{Y,D}_{gt}   — Sieve (ridge on degree-2 basis) predictions
#       fs_pG                 — Pr(G=1|X) logit-lasso
#       fs_{lasso,sieve}_mu_{dgt} — TC conditional means
#       fs_pi_{dgt}           — TC cell propensities
#       ## Cross-fitted predictions (for DR-Wald, DR-TC)
#       cf                    — list with fold_ids, m_{Y,D}_{gt}, pG_raw,
#                               pi_{gt}, mu_Y_{dgt}, pi_{dgt}
#   $lambdas  — pre-registered lambda values (from rep 1 calibration)
#   $n, $M, $K, $p, $timestamp
#
# Usage:
#   Rscript code/R/sim/02_nuisance.R              # all N_VALUES (Lasso, default)
#   Rscript code/R/sim/02_nuisance.R --n 500      # single n
#   Rscript code/R/sim/02_nuisance.R --method rf   # Random Forest
#   Rscript code/R/sim/02_nuisance.R --method nn   # Neural Network
#   Rscript code/R/sim/02_nuisance.R --test       # 20 reps at n=500
#
# =============================================================================

source("code/R/sim/00_config.R")
library(MASS)
library(glmnet)
library(parallel)
library(doParallel)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
TEST_MODE <- "--test" %in% args
n_filter <- NULL
ML_METHOD <- "lasso"  # default
for (i in seq_along(args)) {
  if (args[i] == "--n" && i < length(args)) n_filter <- as.integer(args[i + 1])
  if (args[i] == "--method" && i < length(args)) ML_METHOD <- args[i + 1]
}
n_run <- if (!is.null(n_filter)) n_filter else N_VALUES
M_run <- if (TEST_MODE) 20L else M

# Load method-specific packages
if (ML_METHOD == "rf") { library(randomForest) }
if (ML_METHOD == "nn") { library(nnet) }

# =============================================================================
# ML METHOD DISPATCH: fit a regression or classification model
# =============================================================================
fit_ml_regression <- function(X_train, y_train, X_pred, method = ML_METHOD,
                              lambda_fixed = NULL) {
  n_tr <- nrow(X_train)
  if (method == "lasso") {
    if (!is.null(lambda_fixed)) {
      mod <- glmnet(X_train, y_train, alpha = 1, lambda = lambda_fixed)
    } else {
      mod <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = 5)
    }
    s <- if (!is.null(lambda_fixed)) lambda_fixed else "lambda.min"
    return(as.vector(predict(mod, newx = X_pred, s = s)))
  } else if (method == "rf") {
    mtry <- if (is.null(RF_MTRY)) max(1, floor(ncol(X_train) / 3)) else RF_MTRY
    mod <- randomForest(X_train, y_train, ntree = RF_NTREE, mtry = mtry,
                        nodesize = RF_NODESIZE)
    return(as.vector(predict(mod, newdata = X_pred)))
  } else if (method == "nn") {
    # Standardize inputs for neural network
    mu <- colMeans(X_train); sd_x <- apply(X_train, 2, sd)
    sd_x[sd_x < 1e-8] <- 1
    X_tr_s <- scale(X_train, center = mu, scale = sd_x)
    X_pr_s <- scale(X_pred, center = mu, scale = sd_x)
    mod <- nnet(X_tr_s, y_train, size = NN_SIZE, decay = NN_DECAY,
                maxit = NN_MAXIT, linout = NN_LINOUT, trace = FALSE)
    return(as.vector(predict(mod, newdata = X_pr_s)))
  }
  stop(sprintf("Unknown method: %s", method))
}

fit_ml_classification <- function(X_train, y_train, X_pred, method = ML_METHOD,
                                  lambda_fixed = NULL) {
  # Always use logit-Lasso for propensity scores (RF poorly calibrated)
  if (method %in% c("rf", "nn")) {
    # Use logit-Lasso regardless — RF/NN propensities are unreliable
    if (!is.null(lambda_fixed)) {
      mod <- glmnet(X_train, y_train, alpha = 1, family = "binomial", lambda = lambda_fixed)
    } else {
      mod <- cv.glmnet(X_train, y_train, alpha = 1, family = "binomial", nfolds = 5)
    }
    s <- if (!is.null(lambda_fixed)) lambda_fixed else "lambda.min"
    return(pmax(pmin(as.vector(predict(mod, newx = X_pred, s = s, type = "response")), 0.99), 0.01))
  }
  # Lasso: logit-Lasso
  if (!is.null(lambda_fixed)) {
    mod <- glmnet(X_train, y_train, alpha = 1, family = "binomial", lambda = lambda_fixed)
  } else {
    mod <- cv.glmnet(X_train, y_train, alpha = 1, family = "binomial", nfolds = 5)
  }
  s <- if (!is.null(lambda_fixed)) lambda_fixed else "lambda.min"
  pmax(pmin(as.vector(predict(mod, newx = X_pred, s = s, type = "response")), 0.99), 0.01)
}

# =============================================================================
# PHASE 0: CALIBRATE LAMBDAS (from one representative replication)
# =============================================================================
calibrate_lambdas <- function(W, X_all, X_sieve) {
  n <- nrow(W)
  Y <- W$Y; D <- W$D; G <- W$G; Ti <- W$Ti
  lambdas <- list()

  # m^R_{gt} with Lasso (p=100)
  for (R_name in c("Y", "D")) {
    for (g in 0:1) for (t in 0:1) {
      if (g == 1 && t == 1) next
      key <- sprintf("lasso_%s_%d%d", R_name, g, t)
      idx <- which(G == g & Ti == t)
      lambdas[[key]] <- tryCatch({
        if (length(idx) < 10) 0.01
        else cv.glmnet(X_all[idx, ], W[[R_name]][idx], alpha = 1, nfolds = 5)$lambda.min
      }, error = function(e) 0.01)
    }
  }

  # m^R_{gt} with Sieve (ridge on polynomial basis)
  for (R_name in c("Y", "D")) {
    for (g in 0:1) for (t in 0:1) {
      if (g == 1 && t == 1) next
      key <- sprintf("sieve_%s_%d%d", R_name, g, t)
      idx <- which(G == g & Ti == t)
      lambdas[[key]] <- tryCatch({
        if (length(idx) < 10) 1
        else cv.glmnet(X_sieve[idx, ], W[[R_name]][idx], alpha = 0, nfolds = 5)$lambda.min
      }, error = function(e) 1)
    }
  }

  # Pr(G|X) logit-lasso
  lambdas[["pG"]] <- tryCatch(
    cv.glmnet(X_all, G, alpha = 1, family = "binomial", nfolds = 5)$lambda.min,
    error = function(e) 0.01)

  # TC nuisance: mu^Y_{dgt} with Lasso and Sieve
  for (method in c("lasso", "sieve")) {
    X_use <- if (method == "lasso") X_all else X_sieve
    alpha_val <- if (method == "lasso") 1 else 0
    for (d in 0:1) for (t_val in 0:1) {
      key <- sprintf("tc_%s_%d0%d", method, d, t_val)
      idx <- which(D == d & G == 0 & Ti == t_val)
      lambdas[[key]] <- tryCatch({
        if (length(idx) < 10) ifelse(method == "lasso", 0.01, 1)
        else cv.glmnet(X_use[idx, ], Y[idx], alpha = alpha_val, nfolds = 5)$lambda.min
      }, error = function(e) ifelse(method == "lasso", 0.01, 1))
    }
  }

  # TC propensity: pi_{dgt} logit-lasso
  for (d in 0:1) for (t_val in 0:1) {
    key <- sprintf("pi_tc_%d0%d", d, t_val)
    L_ind <- as.numeric(D == d & G == 0 & Ti == t_val)
    lambdas[[key]] <- tryCatch(
      cv.glmnet(X_all, L_ind, alpha = 1, family = "binomial", nfolds = 5)$lambda.min,
      error = function(e) 0.01)
  }

  lambdas
}

# =============================================================================
# ESTIMATE ALL NUISANCE FOR ONE REPLICATION (fixed lambdas)
# For RF/NN: only cross-fitted nuisance is computed (no full-sample Lasso/Sieve)
# =============================================================================
estimate_nuisance_one <- function(W, X_all, X_sieve, lambdas, K, method = "lasso") {
  n <- nrow(W)
  Y <- W$Y; D <- W$D; G <- W$G; Ti <- W$Ti
  pT <- mean(Ti)
  nuis <- list()

  # ---- Save raw data for Stage 3 ----
  nuis$Y  <- Y
  nuis$D  <- D
  nuis$G  <- G
  nuis$Ti <- Ti
  nuis$V  <- W$V
  nuis$h  <- W$h
  nuis$g0 <- W$g0
  nuis$X  <- X_all

  # ---- Full-sample predictions (Lasso only — Wald_X, TC_X don't use RF/NN) ----
  if (method == "lasso") {
  for (R_name in c("Y", "D")) {
    for (g in 0:1) for (t in 0:1) {
      if (g == 1 && t == 1) next
      key_out <- sprintf("fs_lasso_%s_%d%d", R_name, g, t)
      key_lam <- sprintf("lasso_%s_%d%d", R_name, g, t)
      idx <- which(G == g & Ti == t)
      nuis[[key_out]] <- tryCatch({
        mod <- glmnet(X_all[idx, ], W[[R_name]][idx], alpha = 1, lambda = lambdas[[key_lam]])
        as.vector(predict(mod, newx = X_all))
      }, error = function(e) rep(mean(W[[R_name]][idx]), n))
    }
  }

  # ---- Full-sample Sieve (ridge on degree-2 polynomial basis) ----
  for (R_name in c("Y", "D")) {
    for (g in 0:1) for (t in 0:1) {
      if (g == 1 && t == 1) next
      key_out <- sprintf("fs_sieve_%s_%d%d", R_name, g, t)
      key_lam <- sprintf("sieve_%s_%d%d", R_name, g, t)
      idx <- which(G == g & Ti == t)
      nuis[[key_out]] <- tryCatch({
        mod <- glmnet(X_sieve[idx, ], W[[R_name]][idx], alpha = 0, lambda = lambdas[[key_lam]])
        as.vector(predict(mod, newx = X_sieve))
      }, error = function(e) rep(mean(W[[R_name]][idx]), n))
    }
  }

  # ---- Full-sample Pr(G|X) ----
  nuis$fs_pG <- tryCatch({
    mod <- glmnet(X_all, G, alpha = 1, family = "binomial", lambda = lambdas[["pG"]])
    as.vector(predict(mod, newx = X_all, type = "response"))
  }, error = function(e) rep(mean(G), n))

  # ---- Full-sample TC nuisance (Lasso + Sieve) ----
  for (method in c("lasso", "sieve")) {
    X_use <- if (method == "lasso") X_all else X_sieve
    alpha_val <- if (method == "lasso") 1 else 0
    for (d in 0:1) for (t_val in 0:1) {
      key_out <- sprintf("fs_%s_mu_%d0%d", method, d, t_val)
      key_lam <- sprintf("tc_%s_%d0%d", method, d, t_val)
      idx <- which(D == d & G == 0 & Ti == t_val)
      nuis[[key_out]] <- tryCatch({
        if (length(idx) < 10) rep(if (length(idx) > 0) mean(Y[idx]) else 0, n)
        else {
          mod <- glmnet(X_use[idx, ], Y[idx], alpha = alpha_val, lambda = lambdas[[key_lam]])
          as.vector(predict(mod, newx = X_use))
        }
      }, error = function(e) rep(if (length(idx) > 0) mean(Y[idx]) else 0, n))
    }
  }

  # ---- Full-sample TC propensity ----
  for (d in 0:1) for (t_val in 0:1) {
    key_out <- sprintf("fs_pi_%d0%d", d, t_val)
    key_lam <- sprintf("pi_tc_%d0%d", d, t_val)
    L_ind <- as.numeric(D == d & G == 0 & Ti == t_val)
    nuis[[key_out]] <- tryCatch({
      mod <- glmnet(X_all, L_ind, alpha = 1, family = "binomial", lambda = lambdas[[key_lam]])
      pmax(pmin(as.vector(predict(mod, newx = X_all, type = "response")), 0.99), 0.01)
    }, error = function(e) rep(pmax(mean(L_ind), 0.01), n))
  }
  } # end if (method == "lasso") for full-sample predictions

  # ---- Cross-fitted nuisance (for DR-Wald, DR-TC) — uses method dispatch ----
  fold_ids <- sample(rep(1:K, length.out = n))
  cf <- list(fold_ids = fold_ids)
  cf_names <- c("m_Y_10", "m_Y_01", "m_Y_00", "m_D_10", "m_D_01", "m_D_00",
                "pG_raw", "pi_11", "pi_10", "pi_01", "pi_00",
                "mu_Y_101", "mu_Y_100", "mu_Y_001", "mu_Y_000",
                "pi_101", "pi_100", "pi_001", "pi_000")
  for (nm in cf_names) cf[[nm]] <- numeric(n)

  for (k in 1:K) {
    tr <- which(fold_ids != k); va <- which(fold_ids == k)
    X_tr <- X_all[tr, ]; X_va <- X_all[va, ]

    # Conditional expectations (method-dispatched)
    for (R_name in c("Y", "D")) {
      for (g in 0:1) for (t in 0:1) {
        if (g == 1 && t == 1) next
        key <- paste0("m_", R_name, "_", g, t)
        lam_key <- sprintf("lasso_%s_%d%d", R_name, g, t)
        si <- which(G[tr] == g & Ti[tr] == t)
        if (length(si) < 10) {
          cf[[key]][va] <- mean(W[[R_name]][tr[si]])
          next
        }
        cf[[key]][va] <- tryCatch(
          fit_ml_regression(X_tr[si, , drop = FALSE], W[[R_name]][tr[si]], X_va,
                            method = method, lambda_fixed = lambdas[[lam_key]]),
          error = function(e) rep(mean(W[[R_name]][tr[si]]), length(va)))
      }
    }

    # Propensity Pr(G|X) — always logit-Lasso (via fit_ml_classification)
    cf$pG_raw[va] <- tryCatch(
      fit_ml_classification(X_tr, G[tr], X_va, method = method,
                            lambda_fixed = lambdas[["pG"]]),
      error = function(e) rep(mean(G[tr]), length(va)))

    # Cell probabilities (trimming applied in Stage 3, not here)
    # Save raw pG; Stage 3 applies trimming at any desired level
    pG_va <- cf$pG_raw[va]
    cf$pi_11[va] <- pG_va * pT
    cf$pi_10[va] <- pG_va * (1 - pT)
    cf$pi_01[va] <- (1 - pG_va) * pT
    cf$pi_00[va] <- (1 - pG_va) * (1 - pT)

    # TC nuisance (method-dispatched for mu, logit-Lasso for pi)
    for (d in 0:1) for (t_val in 0:1) {
      key_mu <- paste0("mu_Y_", d, "0", t_val)
      lam_key <- sprintf("tc_lasso_%d0%d", d, t_val)
      si <- which(D[tr] == d & G[tr] == 0 & Ti[tr] == t_val)
      if (length(si) < 10) {
        cf[[key_mu]][va] <- if (length(si) > 0) mean(Y[tr[si]]) else 0
        next
      }
      cf[[key_mu]][va] <- tryCatch(
        fit_ml_regression(X_tr[si, , drop = FALSE], Y[tr[si]], X_va,
                          method = method, lambda_fixed = lambdas[[lam_key]]),
        error = function(e) rep(mean(Y[tr[si]]), length(va)))

      key_pi <- paste0("pi_", d, "0", t_val)
      lam_key_pi <- sprintf("pi_tc_%d0%d", d, t_val)
      L_ind <- as.numeric(D[tr] == d & G[tr] == 0 & Ti[tr] == t_val)
      cf[[key_pi]][va] <- tryCatch(
        fit_ml_classification(X_tr, L_ind, X_va, method = method,
                              lambda_fixed = lambdas[[lam_key_pi]]),
        error = function(e) rep(pmax(mean(L_ind), 0.01), length(va)))
    }
  }

  # Safety bounds on TC propensities
  for (key in c("pi_101", "pi_100", "pi_001", "pi_000")) {
    cf[[key]] <- pmax(pmin(cf[[key]], 0.99), 0.01)
  }

  nuis$cf <- cf
  nuis
}

# =============================================================================
# MAIN LOOP
# =============================================================================
cat("=============================================================================\n")
cat(sprintf("02_nuisance.R %s\n", if (TEST_MODE) "[TEST: M=20]" else ""))
cat(sprintf("  n: %s | M=%d | K=%d | method=%s | cores=%d\n",
            paste(n_run, collapse = ","), M_run, K_FOLDS, ML_METHOD, N_CORES))
cat("=============================================================================\n\n")

for (n_val in n_run) {
  infile <- dgp_file(n_val)
  if (!file.exists(infile)) {
    cat(sprintf("SKIP n=%d: %s not found. Run 01_dgp.R first.\n", n_val, infile))
    next
  }

  cat(sprintf("[%s] === n=%d ===\n", format(Sys.time(), "%H:%M:%S"), n_val))
  raw <- readRDS(infile)
  data_list <- raw$data[1:M_run]

  # Phase 0: calibrate lambdas on rep 1
  cat("  Phase 0: calibrating lambdas...")
  t0 <- Sys.time()
  W1 <- data_list[[1]]
  X_all_1 <- as.matrix(W1[, grep("^X", names(W1))])
  X_sieve_1 <- build_sieve(X_all_1)
  set.seed(42)
  lambdas <- calibrate_lambdas(W1, X_all_1, X_sieve_1)
  cat(sprintf(" done (%.1fs, %d lambdas)\n",
              as.numeric(Sys.time() - t0, units = "secs"), length(lambdas)))

  # Phase 1: estimate nuisance (parallel)
  cat(sprintf("  Phase 1: nuisance estimation (%d reps, %d cores)...\n", M_run, N_CORES))
  t1 <- Sys.time()

  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterExport(cl, c("estimate_nuisance_one", "build_sieve",
                       "fit_ml_regression", "fit_ml_classification",
                       "lambdas", "K_FOLDS", "P_SIEVE", "ML_METHOD",
                       "RF_NTREE", "RF_MTRY", "RF_NODESIZE",
                       "NN_SIZE", "NN_DECAY", "NN_MAXIT", "NN_LINOUT"),
                envir = environment())
  clusterEvalQ(cl, {
    library(glmnet); library(MASS)
    if (ML_METHOD == "rf") library(randomForest)
    if (ML_METHOD == "nn") library(nnet)
  })

  all_nuis <- foreach(r = 1:M_run, .errorhandling = "pass") %dopar% {
    W <- data_list[[r]]
    X_all <- as.matrix(W[, grep("^X", names(W))])
    X_sieve <- build_sieve(X_all)
    estimate_nuisance_one(W, X_all, X_sieve, lambdas, K_FOLDS, method = ML_METHOD)
  }
  stopCluster(cl)

  t1_elapsed <- round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1)
  cat(sprintf("  Phase 1 done: %.1f min\n", t1_elapsed))

  # Count errors
  n_err <- sum(sapply(all_nuis, function(x) inherits(x, "error") || !is.list(x)))
  if (n_err > 0) cat(sprintf("  WARNING: %d/%d reps failed\n", n_err, M_run))

  # Save
  outfile <- nuisance_file(n_val, ML_METHOD)
  saveRDS(list(
    nuisance  = all_nuis,
    lambdas   = lambdas,
    method    = ML_METHOD,
    n = n_val, M = M_run, K = K_FOLDS, p = P_AUG,
    timestamp = Sys.time()
  ), outfile)

  mb <- round(file.size(outfile) / 1e6, 1)
  cat(sprintf("  Saved %s (%.1f MB)\n\n", outfile, mb))
}

cat("=== Stage 2 complete ===\n")
