# Internal helper functions (not exported)

# Minimum training observations for ML fit; below this, return NA
.MIN_CELL_SIZE <- 10L

# Minimum propensity floor — prevents division by zero in IPW weights
.MIN_PROPENSITY <- 0.001

#' Validate inputs for didml estimation
#' @noRd
.validate_inputs <- function(Y, D, G, Ti, X) {
  N <- length(Y)
  if (length(D) != N) stop("`D` must have the same length as `Y`.", call. = FALSE)
  if (length(G) != N) stop("`G` must have the same length as `Y`.", call. = FALSE)
  if (length(Ti) != N) stop("`Ti` must have the same length as `Y`.", call. = FALSE)
  if (nrow(X) != N) stop("`X` must have the same number of rows as `Y`.", call. = FALSE)
  if (!all(G %in% 0:1)) stop("`G` must be binary (0/1).", call. = FALSE)
  if (!all(Ti %in% 0:1)) stop("`Ti` must be binary (0/1).", call. = FALSE)
  if (any(is.na(Y)) || any(is.na(D))) stop("Missing values in `Y` or `D`.", call. = FALSE)
  if (any(is.na(G)) || any(is.na(Ti))) stop("Missing values in `G` or `Ti`.", call. = FALSE)
  if (any(is.na(X))) stop("Missing values in `X`. Remove or impute before estimation.", call. = FALSE)
  invisible(TRUE)
}

#' Create cross-fitting fold assignments (preserves global RNG state)
#' @noRd
.make_folds <- function(N, K, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed))
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else
        rm(".Random.seed", envir = .GlobalEnv)
    })
    set.seed(seed)
  }
  sample(rep(seq_len(K), length.out = N))
}

#' Fit nuisance model within a (g,t) cell using specified method.
#' Returns NA (not cell mean) on failure — never impute.
#' @param family "gaussian" for outcomes, "binomial" for propensity.
#'   Respected by all methods: lasso uses logistic, ols uses glm.fit,
#'   rf clips to [0,1], nn uses logistic output.
#' @noRd
.fit_nuisance_cell <- function(Y_train, X_train, X_predict, method, family = "gaussian") {
  N_pred <- nrow(X_predict)
  na_out <- rep(NA_real_, N_pred)

  # Validate method early
  valid_methods <- c("ols", "lasso", "rf", "nn")
  if (!is.function(method) && !method %in% valid_methods) {
    stop("Unknown method: '", method, "'. Use ",
         paste0("'", valid_methods, "'", collapse = ", "),
         ", or a function.", call. = FALSE)
  }

  if (length(Y_train) < .MIN_CELL_SIZE) {
    warning("Cell has < ", .MIN_CELL_SIZE, " training obs. Returning NA.", call. = FALSE)
    return(na_out)
  }

  tryCatch({
    if (method == "ols") {
      if (family == "binomial") {
        fit <- stats::glm.fit(x = cbind(1, X_train), y = Y_train,
                              family = stats::binomial())
        preds <- as.vector(cbind(1, X_predict) %*% fit$coefficients)
        as.vector(1 / (1 + exp(-preds)))
      } else {
        fit <- stats::lm.fit(x = cbind(1, X_train), y = Y_train)
        coefs <- fit$coefficients
        if (any(is.na(coefs))) {
          warning("Collinear covariates in OLS cell fit. Returning NA.", call. = FALSE)
          return(na_out)
        }
        as.vector(cbind(1, X_predict) %*% coefs)
      }

    } else if (method == "lasso") {
      if (family == "binomial") {
        mod <- glmnet::cv.glmnet(X_train, Y_train, alpha = 1, family = "binomial", nfolds = 5)
        as.vector(stats::predict(mod, newx = X_predict, s = "lambda.min", type = "response"))
      } else {
        mod <- glmnet::cv.glmnet(X_train, Y_train, alpha = 1, nfolds = 5)
        as.vector(stats::predict(mod, newx = X_predict, s = "lambda.min"))
      }

    } else if (method == "rf") {
      if (!requireNamespace("randomForest", quietly = TRUE))
        stop("Install 'randomForest' to use method = 'rf'.", call. = FALSE)
      col_names <- paste0("V", seq_len(ncol(X_train)))
      df_train <- data.frame(y = Y_train, X_train)
      names(df_train)[-1] <- col_names
      mod <- randomForest::randomForest(y ~ ., data = df_train, ntree = 500L,
                                         mtry = max(1L, floor(ncol(X_train) / 3)),
                                         nodesize = 5L)
      df_pred <- data.frame(X_predict)
      names(df_pred) <- col_names
      preds <- as.vector(stats::predict(mod, newdata = df_pred))
      if (family == "binomial") pmax(pmin(preds, 1), 0) else preds

    } else if (method == "nn") {
      if (!requireNamespace("nnet", quietly = TRUE))
        stop("Install 'nnet' to use method = 'nn'.", call. = FALSE)
      X_scaled <- scale(X_train)
      center <- attr(X_scaled, "scaled:center")
      sc <- attr(X_scaled, "scaled:scale")
      sc[sc == 0] <- 1
      is_binary <- (family == "binomial")
      mod <- nnet::nnet(X_scaled, Y_train, size = 5L, decay = 0.01,
                        linout = !is_binary, trace = FALSE, maxit = 500L)
      X_pred_scaled <- scale(X_predict, center = center, scale = sc)
      preds <- as.vector(stats::predict(mod, newdata = X_pred_scaled))
      if (is_binary) pmax(pmin(preds, 1), 0) else preds

    } else if (is.function(method)) {
      method(Y_train, X_train, X_predict)
    }
  }, error = function(e) {
    warning("Nuisance estimation failed in cell: ", e$message, ". Returning NA.", call. = FALSE)
    na_out
  })
}

#' Compute ensemble weights for stacking
#'
#' Given a matrix of cross-validated predictions (one column per learner)
#' and the true outcome, compute combination weights under the specified
#' ensemble strategy.
#'
#' @param Y_cv Numeric vector of true outcomes (length N).
#' @param pred_matrix Numeric matrix of CV predictions (N x J).
#' @param type Character: `"average"`, `"singlebest"`, or `"nnls1"`.
#' @return Numeric vector of weights (length J), summing to 1.
#' @noRd
.compute_ensemble_weights <- function(Y_cv, pred_matrix, type = "average") {
  J <- ncol(pred_matrix)
  if (J == 1L) return(1)

  if (type == "average") return(rep(1 / J, J))

  if (type == "singlebest") {
    mse <- colMeans((Y_cv - pred_matrix)^2, na.rm = TRUE)
    w <- rep(0, J)
    w[which.min(mse)] <- 1
    return(w)
  }

  if (type == "nnls1") {
    # Constrained LS: min ||Y - Xw||^2 s.t. w >= 0, sum(w) = 1
    # Use glmnet ridge with lower.limits = 0, then normalize
    fit <- glmnet::glmnet(pred_matrix, Y_cv, alpha = 0, lambda = 1e-6,
                          lower.limits = 0, intercept = FALSE)
    w <- as.vector(stats::coef(fit))[-1]
    w <- pmax(w, 0)
    if (sum(w) > 0) w <- w / sum(w) else w <- rep(1 / J, J)
    return(w)
  }

  # Fallback
  rep(1 / J, J)
}


#' Fit a stacked (ensemble) nuisance model
#'
#' Calls `.fit_nuisance_cell()` for each method in `methods`, combines
#' predictions using the selected ensemble strategy.
#'
#' @param Y_train Numeric vector of training outcomes.
#' @param X_train Numeric matrix of training covariates.
#' @param X_predict Numeric matrix of prediction covariates.
#' @param methods Character vector of method names (e.g. `c("ols", "lasso")`).
#' @param ensemble_type Character: `"average"`, `"singlebest"`, or `"nnls1"`.
#' @param family Character: `"gaussian"` or `"binomial"`.
#' @param Y_cv Numeric vector of CV outcomes for weight estimation (same as
#'   Y_train for within-fold stacking). Used by `"singlebest"` and `"nnls1"`.
#' @param X_cv Numeric matrix of CV covariates matching Y_cv.
#' @return A list with `predictions` (numeric vector) and `weights` (numeric vector).
#' @noRd
.fit_nuisance_stacked <- function(Y_train, X_train, X_predict, methods,
                                   ensemble_type = "average",
                                   family = "gaussian",
                                   Y_cv = NULL, X_cv = NULL) {
  J <- length(methods)
  N_pred <- nrow(X_predict)

  # Get predictions from each learner
  pred_matrix <- matrix(NA_real_, nrow = N_pred, ncol = J)
  for (j in seq_len(J)) {
    pred_matrix[, j] <- .fit_nuisance_cell(Y_train, X_train, X_predict,
                                            methods[j], family = family)
  }

  # For weight computation: get in-sample (training) predictions per learner
  if (ensemble_type != "average" && !is.null(Y_cv) && !is.null(X_cv)) {
    cv_pred_matrix <- matrix(NA_real_, nrow = length(Y_cv), ncol = J)
    for (j in seq_len(J)) {
      cv_pred_matrix[, j] <- .fit_nuisance_cell(Y_train, X_train, X_cv,
                                                  methods[j], family = family)
    }
    weights <- .compute_ensemble_weights(Y_cv, cv_pred_matrix, type = ensemble_type)
  } else if (ensemble_type != "average") {
    # Fallback: use training predictions on prediction set for weight computation
    weights <- .compute_ensemble_weights(
      rep(NA_real_, N_pred), pred_matrix, type = "average")
  } else {
    weights <- rep(1 / J, J)
  }

  # Combine predictions
  predictions <- as.vector(pred_matrix %*% weights)

  list(predictions = predictions, weights = weights)
}
