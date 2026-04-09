#' Deprecated fddml Functions
#'
#' These functions are deprecated in favor of the \code{didml} family.
#' They will be removed in a future release.
#'
#' @name didml-deprecated
#' @rdname didml-deprecated
#' @param Y Numeric vector of outcomes.
#' @param D Numeric vector of treatment.
#' @param G Binary vector of group assignment (0/1).
#' @param Ti Binary vector of time period (0/1).
#' @param X Numeric matrix of covariates.
#' @param estimand Character: passed as \code{estimator} to \code{didml()}.
#' @param ... Additional arguments passed to the replacement function.
#' @export
fddml <- function(Y, D, G, Ti, X, estimand = "both", ...) {
  .Deprecated("didml", package = "didml",
    msg = "fddml() is deprecated. Use didml(..., iv = TRUE) instead.")
  didml(Y = Y, D = D, G = G, Ti = Ti, X = X,
        design = "2x2", iv = TRUE, dml = TRUE,
        estimator = estimand, ...)
}

#' @rdname didml-deprecated
#' @export
fddml_nuisance <- function(...) {
  .Deprecated("didml_nuisance")
  didml_nuisance(...)
}

#' @rdname didml-deprecated
#' @export
fddml_wald <- function(...) {
  .Deprecated("didml_wald")
  didml_wald(...)
}

#' @rdname didml-deprecated
#' @export
fddml_tc <- function(...) {
  .Deprecated("didml_tc")
  didml_tc(...)
}

#' @rdname didml-deprecated
#' @export
fddml_trim <- function(...) {
  .Deprecated("didml_trim")
  didml_trim(...)
}

#' @rdname didml-deprecated
#' @export
fddml_inference <- function(...) {
  .Deprecated("didml_inference")
  didml_inference(...)
}
