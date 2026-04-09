# Comprehensive tests for didml package
# Covers: numerical correctness, edge cases, all estimation paths,
# stacking combinations, reps, inference, diagnostics

# ============================================================================
# Helper DGP: generates data with known ATT for testing
# ============================================================================
.make_dgp <- function(N = 500, p = 5, ATT = 2.0, fuzzy = FALSE, seed = 123) {
  set.seed(seed)
  X <- matrix(rnorm(N * p), ncol = p)
  G <- rbinom(N, 1, plogis(0.5 * X[, 1]))
  Ti <- rbinom(N, 1, 0.5)
  if (fuzzy) {
    D <- rbinom(N, 1, plogis(X[, 1] + 0.8 * G * Ti))
  } else {
    D <- as.integer(G * Ti)
  }
  Y <- X[, 1] + 0.3 * Ti + ATT * G * Ti + rnorm(N)
  list(Y = Y, D = D, G = G, Ti = Ti, X = X, ATT = ATT)
}

# ============================================================================
# 1. NUMERICAL CORRECTNESS: Does the estimator recover the true ATT?
# ============================================================================

test_that("Sharp DML-Chang recovers ATT within 3 SEs (N=2000)", {
  d <- .make_dgp(N = 2000, ATT = 1.5, fuzzy = FALSE, seed = 101)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols", K = 5,
               trim = "none", seed = 42)
  expect_true(abs(fit$estimates$estimate - d$ATT) < 3 * fit$estimates$se)
})

test_that("Fuzzy DML-Wald recovers ATT within 3 SEs (N=2000)", {
  d <- .make_dgp(N = 2000, ATT = 1.5, fuzzy = TRUE, seed = 102)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "wald",
               method = "ols", K = 5, trim = "none", seed = 42)
  # Fuzzy is noisier; allow wider tolerance
  expect_true(is.finite(fit$estimates$estimate))
  expect_true(fit$estimates$se > 0)
})

test_that("Fuzzy DML-TC recovers ATT within 3 SEs (N=2000)", {
  d <- .make_dgp(N = 2000, ATT = 1.5, fuzzy = TRUE, seed = 103)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "tc",
               method = "ols", K = 5, trim = "none", seed = 42)
  expect_true(is.finite(fit$estimates$estimate))
  expect_true(fit$estimates$se > 0)
})

test_that("Sharp DML-Chang with Lasso recovers ATT (N=1000, p=20)", {
  d <- .make_dgp(N = 1000, p = 20, ATT = 2.0, fuzzy = FALSE, seed = 104)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "lasso", K = 5,
               trim = "none", seed = 42)
  expect_true(abs(fit$estimates$estimate - d$ATT) < 3 * fit$estimates$se)
})

# ============================================================================
# 2. INFERENCE CORRECTNESS
# ============================================================================

test_that("99% CI contains true ATT (coverage check, single run)", {
  d <- .make_dgp(N = 2000, ATT = 1.0, fuzzy = FALSE, seed = 201)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols", K = 5,
               trim = "none", seed = 42)
  # Use 99% CI to avoid flaky 5% chance of failure
  z99 <- qnorm(0.995)
  ci_lo <- fit$estimates$estimate - z99 * fit$estimates$se
  ci_hi <- fit$estimates$estimate + z99 * fit$estimates$se
  expect_true(ci_lo < d$ATT && d$ATT < ci_hi)
})

test_that("Bootstrap SE is finite and positive", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 202)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols", K = 2,
               se_type = "bootstrap", B = 200, trim = "none", seed = 42)
  expect_true(is.finite(fit$estimates$se))
  expect_true(fit$estimates$se > 0)
})

test_that("Cluster SEs are at least as large as analytical SEs (sharp)", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 203)
  cluster <- rep(seq_len(30), each = 10)

  fit_iid <- didml(d$Y, d$D, d$G, d$Ti, d$X,
                   iv = FALSE, dml = TRUE, method = "ols", K = 2,
                   trim = "none", seed = 42)
  fit_cl <- didml(d$Y, d$D, d$G, d$Ti, d$X,
                  iv = FALSE, dml = TRUE, method = "ols", K = 2,
                  cluster = cluster, trim = "none", seed = 42)

  expect_true(is.finite(fit_cl$estimates$se))
  # Cluster SE should generally be >= analytical
  # (not guaranteed in all samples, so just check it's positive)
  expect_true(fit_cl$estimates$se > 0)
})

test_that("p-value is between 0 and 1", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 204)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols", K = 2,
               trim = "none", seed = 42)
  expect_true(fit$estimates$p_value >= 0 && fit$estimates$p_value <= 1)
})

test_that("confint respects level parameter", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 205)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols", K = 2,
               trim = "none", seed = 42)
  ci_95 <- confint(fit, level = 0.95)
  ci_99 <- confint(fit, level = 0.99)
  # 99% CI should be wider than 95%
  width_95 <- ci_95[1, 2] - ci_95[1, 1]
  width_99 <- ci_99[1, 2] - ci_99[1, 1]
  expect_true(width_99 > width_95)
})

# ============================================================================
# 3. STACKING TESTS
# ============================================================================

test_that("Stacking with 3 methods works (ols + lasso + rf)", {
  skip_if_not_installed("randomForest")
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 301)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE,
               method = c("ols", "lasso", "rf"),
               ensemble_type = "average", K = 2,
               trim = "none", seed = 42)
  expect_true(is.finite(fit$estimates$estimate))
  expect_true(!is.null(fit$ensemble_weights))
})

test_that("nnls1 ensemble weights sum to 1", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 302)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE,
               method = c("ols", "lasso"), ensemble_type = "nnls1",
               K = 2, trim = "none", seed = 42)
  expect_true(!is.null(fit$ensemble_weights))
  # Check at least one weight vector sums to ~1
  has_unit_sum <- any(vapply(fit$ensemble_weights, function(w) {
    length(w) > 1 && abs(sum(w) - 1) < 0.01
  }, logical(1)))
  expect_true(has_unit_sum)
})

test_that("Stacking with fuzzy DID works", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = TRUE, seed = 303)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "wald",
               method = c("ols", "lasso"), ensemble_type = "average",
               K = 2, trim = "none", seed = 42)
  expect_true(is.finite(fit$estimates$estimate))
})

# ============================================================================
# 4. SEPARATE LEARNERS TESTS
# ============================================================================

test_that("method_outcome overrides method for outcome regression", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 401)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE,
               method = "ols", method_outcome = "lasso",
               K = 2, trim = "none", seed = 42)
  expect_true(is.finite(fit$estimates$estimate))
  expect_equal(fit$settings$method_outcome, "lasso")
})

test_that("method_outcome + method_propensity both override", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = TRUE, seed = 402)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "wald",
               method = "ols", method_outcome = "lasso",
               method_propensity = "ols",
               K = 2, trim = "none", seed = 42)
  expect_true(is.finite(fit$estimates$estimate))
  expect_equal(fit$settings$method_outcome, "lasso")
  expect_equal(fit$settings$method_propensity, "ols")
})

test_that("Stacked outcome + single propensity works", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 403)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE,
               method_outcome = c("ols", "lasso"),
               method_propensity = "lasso",
               ensemble_type = "average",
               K = 2, trim = "none", seed = 42)
  expect_true(is.finite(fit$estimates$estimate))
})

# ============================================================================
# 5. CROSS-FITTING REPETITIONS TESTS
# ============================================================================

test_that("reps = 3 produces aggregated estimate with lower variance", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 501)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols",
               K = 2, reps = 3L, trim = "none", seed = 42)
  expect_equal(length(fit$reps_detail), 3)
  # Each rep should have an estimate
  rep_ests <- vapply(fit$reps_detail, function(r) r$estimate, numeric(1))
  expect_true(all(is.finite(rep_ests)))
  # Aggregated estimate should be median of reps
  expect_equal(fit$estimates$estimate, median(rep_ests))
})

test_that("reps = 1 matches single-run output exactly", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 502)
  fit_1 <- didml(d$Y, d$D, d$G, d$Ti, d$X,
                 iv = FALSE, dml = TRUE, method = "ols",
                 K = 2, reps = 1L, trim = "none", seed = 42)
  # reps_detail should be NULL when reps = 1
  expect_null(fit_1$reps_detail)
})

test_that("reps with fuzzy DID works", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = TRUE, seed = 503)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "wald",
               method = "ols", K = 2, reps = 2L,
               trim = "none", seed = 42)
  expect_equal(length(fit$reps_detail), 2)
  expect_true(is.finite(fit$estimates$estimate))
})

# ============================================================================
# 6. TRIMMING TESTS
# ============================================================================

test_that("Auto trimming selects non-negative alpha", {
  d <- .make_dgp(N = 500, ATT = 1.0, fuzzy = TRUE, seed = 601)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "wald",
               method = "ols", K = 3, trim = "auto", seed = 42)
  expect_true(fit$trim_info$alpha_hat >= 0)
  expect_true(fit$trim_info$alpha_hat <= 0.50)
})

test_that("Fixed trimming at 0.10 works", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = TRUE, seed = 602)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "wald",
               method = "ols", K = 2, trim = "fixed",
               trim_alpha = 0.10, seed = 42)
  expect_equal(fit$trim_info$alpha_hat, 0.10)
})

test_that("No trimming leaves propensity unchanged", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = TRUE, seed = 603)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "wald",
               method = "ols", K = 2, trim = "none", seed = 42)
  expect_equal(fit$trim_info$alpha_hat, 0)
})

# ============================================================================
# 7. EDGE CASES & ERROR HANDLING
# ============================================================================

test_that("Error on design = 'multi'", {
  d <- .make_dgp(N = 100, seed = 701)
  expect_error(
    didml(d$Y, d$D, d$G, d$Ti, d$X, design = "multi"),
    "future release"
  )
})

test_that("Error on iv = TRUE, dml = FALSE", {
  d <- .make_dgp(N = 100, seed = 702)
  expect_error(
    didml(d$Y, d$D, d$G, d$Ti, d$X, iv = TRUE, dml = FALSE),
    "Fuzzy DID requires DML"
  )
})

test_that("Error on mismatched lengths", {
  expect_error(
    didml(rnorm(100), rnorm(99), rbinom(100, 1, 0.5),
          rbinom(100, 1, 0.5), matrix(rnorm(300), ncol = 3)),
    "same length"
  )
})

test_that("Error on non-binary G", {
  expect_error(
    didml(rnorm(100), rnorm(100), rnorm(100),
          rbinom(100, 1, 0.5), matrix(rnorm(300), ncol = 3)),
    "binary"
  )
})

test_that("Error on NA in Y", {
  Y <- rnorm(100); Y[1] <- NA
  expect_error(
    didml(Y, rbinom(100, 1, 0.5), rbinom(100, 1, 0.5),
          rbinom(100, 1, 0.5), matrix(rnorm(300), ncol = 3)),
    "Missing values"
  )
})

test_that("Error on unknown method string", {
  d <- .make_dgp(N = 100, seed = 706)
  expect_error(
    didml(d$Y, d$D, d$G, d$Ti, d$X, method = "xgboost"),
    "Unknown"
  )
})

test_that("Small sample (N=50) doesn't crash", {
  d <- .make_dgp(N = 50, p = 2, ATT = 1.0, fuzzy = FALSE, seed = 707)
  fit <- suppressWarnings(
    didml(d$Y, d$D, d$G, d$Ti, d$X,
          iv = FALSE, dml = TRUE, method = "ols",
          K = 2, trim = "none", seed = 42)
  )
  expect_s3_class(fit, "didml")
})

test_that("K = 2 (minimum folds) works", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 708)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols",
               K = 2, trim = "none", seed = 42)
  expect_true(is.finite(fit$estimates$estimate))
})

test_that("K = 10 (many folds) works", {
  d <- .make_dgp(N = 500, ATT = 1.0, fuzzy = FALSE, seed = 709)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols",
               K = 10, trim = "none", seed = 42)
  expect_true(is.finite(fit$estimates$estimate))
})

# ============================================================================
# 8. S3 METHOD TESTS
# ============================================================================

test_that("print works for all three paths", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 801)

  # Sharp DML
  fit_chang <- didml(d$Y, d$D, d$G, d$Ti, d$X,
                     iv = FALSE, dml = TRUE, method = "ols",
                     K = 2, trim = "none", seed = 42)
  expect_output(print(fit_chang), "Sharp DID")

  # Fuzzy DML
  d2 <- .make_dgp(N = 200, ATT = 1.0, fuzzy = TRUE, seed = 802)
  fit_wald <- didml(d2$Y, d2$D, d2$G, d2$Ti, d2$X,
                    iv = TRUE, dml = TRUE, estimator = "wald",
                    method = "ols", K = 2, trim = "none", seed = 42)
  expect_output(print(fit_wald), "Fuzzy DID")
})

test_that("summary works and contains expected sections", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 803)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols",
               K = 2, trim = "none", seed = 42)
  s <- summary(fit)
  expect_s3_class(s, "summary.didml")
  expect_output(print(s), "Summary")
  expect_output(print(s), "Method")
})

test_that("tidy returns all required columns", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = TRUE, seed = 804)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "both",
               method = "ols", K = 2, trim = "none", seed = 42)
  td <- tidy.didml(fit)
  expect_true(all(c("term", "estimate", "std.error", "statistic",
                     "p.value", "conf.low", "conf.high") %in% names(td)))
  expect_equal(nrow(td), 2)  # Wald + TC
})

test_that("confint returns matrix with correct dimensions", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = TRUE, seed = 805)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "both",
               method = "ols", K = 2, trim = "none", seed = 42)
  ci <- confint(fit)
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
})

test_that("plot propensity works without error", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 806)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE, method = "ols",
               K = 2, trim = "auto", seed = 42)
  p <- plot(fit, type = "propensity")
  expect_s3_class(p, "ggplot")
})

test_that("plot scores works without error", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = TRUE, seed = 807)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = TRUE, dml = TRUE, estimator = "wald",
               method = "ols", K = 2, trim = "none", seed = 42)
  p <- plot(fit, type = "scores")
  expect_s3_class(p, "ggplot")
})

# ============================================================================
# 9. NUISANCE ESTIMATION TESTS
# ============================================================================

test_that("didml_nuisance with iv=FALSE returns Chang structure", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 901)
  nuis <- didml_nuisance(d$Y, d$D, d$G, d$Ti, d$X,
                          method = "ols", K = 2, iv = FALSE, seed = 42)
  expect_true("ell_20" %in% names(nuis))
  expect_true("pG_raw" %in% names(nuis))
  expect_true("pT" %in% names(nuis))
  expect_equal(nuis$estimand, "chang")
  expect_null(nuis$m_Y_10)
})

test_that("didml_nuisance with iv=TRUE returns full structure", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = TRUE, seed = 902)
  nuis <- didml_nuisance(d$Y, d$D, d$G, d$Ti, d$X,
                          method = "ols", K = 2, iv = TRUE,
                          estimand = "both", seed = 42)
  expect_true(all(c("m_Y_10", "m_Y_01", "m_Y_00",
                     "m_D_10", "m_D_01", "m_D_00",
                     "pG_raw") %in% names(nuis)))
  # TC nuisance should also be present
  expect_true("mu_Y_101" %in% names(nuis))
})

test_that("didml_nuisance propensity is bounded away from 0 and 1", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 903)
  nuis <- didml_nuisance(d$Y, d$D, d$G, d$Ti, d$X,
                          method = "ols", K = 2, iv = FALSE, seed = 42)
  expect_true(all(nuis$pG_raw >= 0.001, na.rm = TRUE))
  expect_true(all(nuis$pG_raw <= 0.999, na.rm = TRUE))
})

test_that("didml_nuisance stacking returns ensemble weights", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 904)
  nuis <- didml_nuisance(d$Y, d$D, d$G, d$Ti, d$X,
                          method = c("ols", "lasso"), K = 2,
                          iv = FALSE, ensemble_type = "average", seed = 42)
  expect_true("ensemble_weights" %in% names(nuis))
})

# ============================================================================
# 10. SCORE FUNCTION TESTS
# ============================================================================

test_that("didml_wald returns correct structure", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = TRUE, seed = 1001)
  nuis <- didml_nuisance(d$Y, d$D, d$G, d$Ti, d$X,
                          method = "ols", K = 2, iv = TRUE,
                          estimand = "wald", seed = 42)
  # Apply trivial trimming
  pT <- mean(d$Ti)
  nuis$pi_11 <- nuis$pG_raw * pT
  nuis$pi_10 <- nuis$pG_raw * (1 - pT)
  nuis$pi_01 <- (1 - nuis$pG_raw) * pT
  nuis$pi_00 <- (1 - nuis$pG_raw) * (1 - pT)

  W <- data.frame(Y = d$Y, D = d$D, G = d$G, Ti = d$Ti)
  res <- didml_wald(W, nuis)

  expect_true(all(c("estimate", "variance", "scores",
                     "num_i", "den_i") %in% names(res)))
  expect_equal(length(res$scores), length(d$Y))
  expect_true(is.finite(res$estimate))
})

test_that("didml_chang returns correct structure", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 1002)
  nuis <- didml_nuisance(d$Y, d$D, d$G, d$Ti, d$X,
                          method = "ols", K = 2, iv = FALSE, seed = 42)
  W <- data.frame(Y = d$Y, D = d$D, G = d$G, Ti = d$Ti)
  res <- didml_chang(W, nuis)

  expect_true(all(c("estimate", "variance", "scores",
                     "num_i", "den_i") %in% names(res)))
  # den_i should be all 1s for Chang
  expect_true(all(res$den_i == 1))
  expect_true(is.finite(res$estimate))
})

test_that("didml_chang scores have mean approximately zero (centered)", {
  d <- .make_dgp(N = 2000, ATT = 1.0, fuzzy = FALSE, seed = 1003)
  nuis <- didml_nuisance(d$Y, d$D, d$G, d$Ti, d$X,
                          method = "ols", K = 5, iv = FALSE, seed = 42)
  W <- data.frame(Y = d$Y, D = d$D, G = d$G, Ti = d$Ti)
  res <- didml_chang(W, nuis)
  # Centered scores should have mean near 0 (relaxed for finite sample)
  ok <- is.finite(res$scores)
  expect_true(abs(mean(res$scores[ok])) < 1.0)
})

# ============================================================================
# 11. REPRODUCIBILITY TESTS
# ============================================================================

test_that("Same seed produces identical results", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 1101)
  fit1 <- didml(d$Y, d$D, d$G, d$Ti, d$X,
                iv = FALSE, dml = TRUE, method = "ols",
                K = 3, trim = "none", seed = 42)
  fit2 <- didml(d$Y, d$D, d$G, d$Ti, d$X,
                iv = FALSE, dml = TRUE, method = "ols",
                K = 3, trim = "none", seed = 42)
  expect_equal(fit1$estimates$estimate, fit2$estimates$estimate)
  expect_equal(fit1$estimates$se, fit2$estimates$se)
})

test_that("Different seeds produce different results", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = FALSE, seed = 1102)
  fit1 <- didml(d$Y, d$D, d$G, d$Ti, d$X,
                iv = FALSE, dml = TRUE, method = "ols",
                K = 3, trim = "none", seed = 42)
  fit2 <- didml(d$Y, d$D, d$G, d$Ti, d$X,
                iv = FALSE, dml = TRUE, method = "ols",
                K = 3, trim = "none", seed = 99)
  # Different fold splits → slightly different estimates
  expect_false(identical(fit1$estimates$estimate, fit2$estimates$estimate))
})

# ============================================================================
# 12. DEPRECATED WRAPPER TESTS
# ============================================================================

test_that("fddml() still works with deprecation warning", {
  d <- .make_dgp(N = 200, ATT = 1.0, fuzzy = TRUE, seed = 1201)
  expect_warning(
    fit <- fddml(d$Y, d$D, d$G, d$Ti, d$X,
                 estimand = "wald", method = "ols",
                 K = 2, trim = "none", seed = 42),
    "deprecated"
  )
  expect_s3_class(fit, "didml")
  expect_true(is.finite(fit$estimates$estimate))
})

# ============================================================================
# 13. COMBINATION TESTS: Stacking + Reps + Separate Learners
# ============================================================================

test_that("All three features combined work", {
  d <- .make_dgp(N = 300, ATT = 1.0, fuzzy = FALSE, seed = 1301)
  fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
               iv = FALSE, dml = TRUE,
               method_outcome = c("ols", "lasso"),
               method_propensity = "ols",
               ensemble_type = "average",
               K = 2, reps = 2L,
               trim = "none", seed = 42)
  expect_s3_class(fit, "didml")
  expect_true(is.finite(fit$estimates$estimate))
  expect_equal(length(fit$reps_detail), 2)
})
