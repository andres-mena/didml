# Tests for didml package — new dispatcher and Chang estimator

test_that("didml sharp DML-DID (Chang) returns correct structure", {
  set.seed(42)
  N <- 200
  X <- matrix(rnorm(N * 3), ncol = 3)
  G <- rbinom(N, 1, plogis(0.3 * X[, 1]))
  Ti <- rbinom(N, 1, 0.5)
  D <- as.integer(G * Ti)  # sharp: perfect compliance
  Y <- 1 + X[, 1] + 0.5 * G * Ti + rnorm(N)

  fit <- didml(Y, D, G, Ti, X, iv = FALSE, dml = TRUE,
               method = "ols", K = 2, trim = "none", seed = 42)

  expect_s3_class(fit, "didml")
  expect_s3_class(fit, "didml_sharp")
  expect_equal(nrow(fit$estimates), 1)
  expect_true(grepl("Chang", fit$estimates$estimand))
  expect_true(is.finite(fit$estimates$estimate))
  expect_true(fit$estimates$se > 0)
})

test_that("didml sharp DML-DID recovers known ATT", {
  set.seed(123)
  N <- 2000
  ATT_true <- 2.0
  X <- matrix(rnorm(N * 5), ncol = 5)
  G <- rbinom(N, 1, plogis(0.5 * X[, 1]))
  Ti <- rbinom(N, 1, 0.5)
  D <- as.integer(G * Ti)
  Y <- X[, 1] + 0.3 * Ti + ATT_true * G * Ti + rnorm(N)

  fit <- didml(Y, D, G, Ti, X, iv = FALSE, dml = TRUE,
               method = "ols", K = 5, trim = "none", seed = 42)

  expect_true(abs(fit$estimates$estimate - ATT_true) < 3 * fit$estimates$se)
})

test_that("didml fuzzy Wald returns correct structure", {
  set.seed(42)
  N <- 200
  X <- matrix(rnorm(N * 3), ncol = 3)
  G <- rbinom(N, 1, plogis(0.3 * X[, 1]))
  Ti <- rbinom(N, 1, 0.5)
  D <- rbinom(N, 1, plogis(X[, 1] + 0.5 * G * Ti))
  Y <- 1 + X[, 1] + 0.3 * Ti + D * G * Ti + rnorm(N)

  fit <- didml(Y, D, G, Ti, X, iv = TRUE, dml = TRUE,
               estimator = "wald", method = "ols",
               K = 2, trim = "none", seed = 42)

  expect_s3_class(fit, "didml")
  expect_s3_class(fit, "didml_fuzzy")
  expect_equal(nrow(fit$estimates), 1)
  expect_equal(fit$estimates$estimand, "DML-Wald")
  expect_true(is.finite(fit$estimates$estimate))
})

test_that("didml fuzzy both returns two rows", {
  set.seed(42)
  N <- 200
  X <- matrix(rnorm(N * 3), ncol = 3)
  G <- rbinom(N, 1, 0.5)
  Ti <- rbinom(N, 1, 0.5)
  D <- rbinom(N, 1, plogis(X[, 1] + G))
  Y <- 1 + X[, 1] + D + rnorm(N)

  fit <- didml(Y, D, G, Ti, X, iv = TRUE, dml = TRUE,
               estimator = "both", method = "ols",
               K = 2, trim = "none", seed = 42)

  expect_equal(nrow(fit$estimates), 2)
  expect_true("DML-Wald" %in% fit$estimates$estimand)
  expect_true("DML-TC" %in% fit$estimates$estimand)
})

test_that("didml multi-period errors with informative message", {
  expect_error(
    didml(rnorm(10), rbinom(10, 1, 0.5), rbinom(10, 1, 0.5),
          rbinom(10, 1, 0.5), matrix(rnorm(30), ncol = 3),
          design = "multi"),
    "future release"
  )
})

test_that("didml fuzzy without DML errors", {
  expect_error(
    didml(rnorm(10), rbinom(10, 1, 0.5), rbinom(10, 1, 0.5),
          rbinom(10, 1, 0.5), matrix(rnorm(30), ncol = 3),
          iv = TRUE, dml = FALSE),
    "Fuzzy DID requires DML"
  )
})

test_that("didml Chang nuisance has correct structure", {
  set.seed(42)
  N <- 200
  X <- matrix(rnorm(N * 3), ncol = 3)
  G <- rbinom(N, 1, 0.5)
  Ti <- rbinom(N, 1, 0.5)
  Y <- 1 + X[, 1] + rnorm(N)
  D <- as.integer(G * Ti)

  nuis <- didml_nuisance(Y, D, G, Ti, X, method = "ols", K = 2,
                          iv = FALSE, seed = 42)

  expect_true("ell_20" %in% names(nuis))
  expect_true("pG_raw" %in% names(nuis))
  expect_true("pT" %in% names(nuis))
  expect_equal(nuis$estimand, "chang")
  # Should NOT have m_Y or m_D fields
  expect_null(nuis$m_Y_10)
  expect_null(nuis$m_D_10)
})

test_that("didml print and summary methods work", {
  set.seed(42)
  N <- 100
  X <- matrix(rnorm(N * 2), ncol = 2)
  G <- rbinom(N, 1, 0.5)
  Ti <- rbinom(N, 1, 0.5)
  D <- as.integer(G * Ti)
  Y <- rnorm(N)

  fit <- didml(Y, D, G, Ti, X, iv = FALSE, dml = TRUE,
               method = "ols", K = 2, trim = "none", seed = 1)

  expect_output(print(fit), "Sharp DID")
  expect_output(print(summary(fit)), "Summary")
})

test_that("didml cluster SEs work for Chang estimator", {
  set.seed(42)
  N <- 300
  n_clust <- 30
  cluster <- rep(seq_len(n_clust), each = N / n_clust)
  X <- matrix(rnorm(N * 3), ncol = 3)
  G <- rbinom(N, 1, 0.5)
  Ti <- rbinom(N, 1, 0.5)
  D <- as.integer(G * Ti)
  Y <- 1 + X[, 1] + D + rnorm(N)

  fit <- didml(Y, D, G, Ti, X, iv = FALSE, dml = TRUE,
               method = "ols", K = 2, trim = "none",
               cluster = cluster, seed = 42)

  expect_true(is.finite(fit$estimates$se))
})
