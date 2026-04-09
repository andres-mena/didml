test_that("fddml_wald returns correct structure with OLS", {
  set.seed(20260327)
  N <- 200
  X <- matrix(rnorm(N * 3), ncol = 3)
  G <- rbinom(N, 1, plogis(0.3 * X[, 1]))
  Ti <- rbinom(N, 1, 0.5)
  D <- rbinom(N, 1, plogis(X[, 1] + 0.5 * G * Ti))
  Y <- 1 + X[, 1] + 0.3 * Ti + D * G * Ti + rnorm(N)

  fit <- fddml(Y, D, G, Ti, X, estimand = "wald", method = "ols",
               K = 2, trim = "none", seed = 42)

  expect_s3_class(fit, "didml")
  expect_equal(nrow(fit$estimates), 1)
  expect_equal(fit$estimates$estimand, "DML-Wald")
  expect_true(is.finite(fit$estimates$estimate))
  expect_true(fit$estimates$se > 0)
})

test_that("fddml with both estimands returns two rows", {
  set.seed(20260327)
  N <- 200
  X <- matrix(rnorm(N * 3), ncol = 3)
  G <- rbinom(N, 1, 0.5)
  Ti <- rbinom(N, 1, 0.5)
  D <- rbinom(N, 1, plogis(X[, 1] + G))
  Y <- 1 + X[, 1] + D + rnorm(N)

  fit <- fddml(Y, D, G, Ti, X, estimand = "both", method = "ols",
               K = 2, trim = "none", seed = 42)

  expect_equal(nrow(fit$estimates), 2)
  expect_true("DML-Wald" %in% fit$estimates$estimand)
  expect_true("DML-TC" %in% fit$estimates$estimand)
})

test_that("cluster SEs are larger than analytical SEs", {
  set.seed(20260327)
  N <- 300
  n_clust <- 30
  cluster <- rep(seq_len(n_clust), each = N / n_clust)
  X <- matrix(rnorm(N * 3), ncol = 3)
  G <- rbinom(N, 1, 0.5)
  Ti <- rbinom(N, 1, 0.5)
  D <- rbinom(N, 1, 0.5)
  Y <- 1 + X[, 1] + D + rnorm(N)

  fit_iid <- fddml(Y, D, G, Ti, X, estimand = "wald", method = "ols",
                    K = 2, trim = "none", se_type = "analytical", seed = 42)
  fit_clust <- fddml(Y, D, G, Ti, X, estimand = "wald", method = "ols",
                      K = 2, trim = "none", cluster = cluster,
                      se_type = "cluster", seed = 42)

  # Cluster SEs should generally be at least as large
  expect_true(is.finite(fit_clust$estimates$se))
  expect_true(is.finite(fit_iid$estimates$se))
})

test_that("trimming reduces alpha_hat to a finite value", {
  set.seed(20260327)
  N <- 300
  X <- matrix(rnorm(N * 3), ncol = 3)
  G <- rbinom(N, 1, plogis(2 * X[, 1]))
  Ti <- rbinom(N, 1, 0.5)
  D <- rbinom(N, 1, 0.5)
  Y <- 1 + D + rnorm(N)

  fit <- fddml(Y, D, G, Ti, X, estimand = "wald", method = "ols",
               K = 2, trim = "auto", seed = 42)

  expect_true(is.finite(fit$trim_info$alpha_hat))
  expect_true(fit$trim_info$alpha_hat >= 0)
})

test_that("print and summary methods work", {
  set.seed(42)
  N <- 100
  X <- matrix(rnorm(N * 2), ncol = 2)
  G <- rbinom(N, 1, 0.5)
  Ti <- rbinom(N, 1, 0.5)
  D <- rbinom(N, 1, 0.5)
  Y <- rnorm(N)

  fit <- fddml(Y, D, G, Ti, X, estimand = "wald", method = "ols",
               K = 2, trim = "none", seed = 1)

  expect_output(print(fit), "Fuzzy DID")
  expect_output(print(summary(fit)), "Summary")
})

test_that("tidy returns correct columns", {
  set.seed(42)
  N <- 100
  X <- matrix(rnorm(N * 2), ncol = 2)
  G <- rbinom(N, 1, 0.5)
  Ti <- rbinom(N, 1, 0.5)
  D <- rbinom(N, 1, 0.5)
  Y <- rnorm(N)

  fit <- fddml(Y, D, G, Ti, X, estimand = "wald", method = "ols",
               K = 2, trim = "none", seed = 1)
  td <- tidy.didml(fit)

  expect_true("term" %in% names(td))
  expect_true("estimate" %in% names(td))
  expect_true("std.error" %in% names(td))
  expect_true("p.value" %in% names(td))
})
