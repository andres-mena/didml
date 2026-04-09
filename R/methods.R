#' Print a didml object
#'
#' @param x A \code{didml} object.
#' @param digits Number of significant digits (default 4).
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @method print didml
#' @export
print.didml <- function(x, digits = 4, ...) {
  # Choose header based on design type

  if (!x$iv && !x$dml) {
    header <- "Sharp DID (DRDID)"
  } else if (!x$iv && x$dml) {
    header <- "Sharp DID with Machine Learning (DML-DID)"
  } else {
    header <- "Fuzzy DID with Machine Learning"
  }

  cat("\n  ", header, "\n\n", sep = "")
  cat("  N =", x$N, " p =", x$p)
  if (x$dml) {
    cat("  K =", x$settings$K, " method =", x$settings$method)
  }
  cat("\n")
  if (!is.null(x$trim_info) && x$trim_info$alpha_hat > 0)
    cat("  Trimming: alpha =", round(x$trim_info$alpha_hat, 3), "\n")
  cat("\n")
  est <- x$estimates
  for (i in seq_len(nrow(est))) {
    if (is.finite(est$estimate[i])) {
      cat(sprintf("  %-10s  %s (SE = %s)  95%% CI [%s, %s]\n",
                  est$estimand[i],
                  formatC(est$estimate[i], digits = digits, format = "f"),
                  formatC(est$se[i], digits = digits, format = "f"),
                  formatC(est$ci_lower[i], digits = digits, format = "f"),
                  formatC(est$ci_upper[i], digits = digits, format = "f")))
    } else {
      cat(sprintf("  %-10s  NA (estimation failed)\n", est$estimand[i]))
    }
  }
  cat("\n")
  invisible(x)
}

#' Summary of a didml object
#'
#' @param object A \code{didml} object.
#' @param ... Additional arguments (ignored).
#' @return An object of class \code{summary.didml}, printed via
#'   \code{\link{print.summary.didml}}.
#' @method summary didml
#' @export
summary.didml <- function(object, ...) {
  out <- list(
    call = object$call,
    N = object$N,
    p = object$p,
    estimates = object$estimates,
    settings = object$settings,
    trim_info = object$trim_info,
    iv = object$iv,
    dml = object$dml,
    wald_naive = if (!is.null(object$wald) && !is.null(object$wald$wald_naive))
      object$wald$wald_naive else NULL
  )
  class(out) <- "summary.didml"
  out
}

#' Print a summary.didml object
#'
#' @param x A \code{summary.didml} object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @method print summary.didml
#' @export
print.summary.didml <- function(x, ...) {
  # Choose header
  if (!x$iv && !x$dml) {
    header <- "Sharp DID (DRDID) - Summary"
  } else if (!x$iv && x$dml) {
    header <- "Sharp DID with Machine Learning (DML-DID) - Summary"
  } else {
    header <- "Fuzzy DID with Machine Learning - Summary"
  }

  cat("\n  ", header, "\n", sep = "")
  cat("  ", paste(rep("-", 50), collapse = ""), "\n")
  cat("  Call:", deparse(x$call, width.cutoff = 60), "\n")
  cat("\n  Sample: N =", x$N, ", p =", x$p, "\n")
  if (x$dml) {
    cat("  Method:", x$settings$method, "| Folds:", x$settings$K, "\n")
  }
  cat("  SE type:", x$settings$se_type, "\n")
  if (!is.null(x$trim_info)) {
    cat("  Trimming:", x$settings$trim,
        if (x$trim_info$alpha_hat > 0)
          paste0("(alpha = ", round(x$trim_info$alpha_hat, 3), ")")
        else "",
        "\n")
  }
  cat("\n  Estimates:\n")
  print(x$estimates, row.names = FALSE, digits = 4)
  cat("\n")
  if (!is.null(x$wald_naive) && is.finite(x$wald_naive))
    cat("  Wald naive (no DML):", round(x$wald_naive, 4), "\n")
  cat("\n")
  invisible(x)
}

#' @method confint didml
#' @export
confint.didml <- function(object, parm = NULL, level = 0.95, ...) {
  z <- stats::qnorm((1 + level) / 2)
  est <- object$estimates
  if (!is.null(parm)) {
    est <- est[est$estimand %in% parm, , drop = FALSE]
    if (nrow(est) == 0)
      stop("No matching estimand for parm = '",
           paste(parm, collapse = "', '"), "'.",
           call. = FALSE)
  }
  ci <- data.frame(
    lower = est$estimate - z * est$se,
    upper = est$estimate + z * est$se,
    row.names = est$estimand
  )
  names(ci) <- paste0(
    formatC(100 * c((1 - level) / 2, (1 + level) / 2),
            format = "f", digits = 1), " %")
  ci
}

#' Tidy a didml object (broom-style)
#'
#' @param x A \code{didml} object.
#' @param ... Additional arguments (ignored).
#' @return A data frame with columns: term, estimate, std.error, statistic,
#'   p.value, conf.low, conf.high.
#' @method tidy didml
#' @export
tidy.didml <- function(x, ...) {
  est <- x$estimates
  data.frame(
    term = est$estimand,
    estimate = est$estimate,
    std.error = est$se,
    statistic = est$estimate / est$se,
    p.value = est$p_value,
    conf.low = est$ci_lower,
    conf.high = est$ci_upper,
    stringsAsFactors = FALSE
  )
}

#' Plot Diagnostics for didml Results
#'
#' Produces diagnostic plots for DID-ML estimation results.
#'
#' @param x A \code{didml} object.
#' @param type Character: \code{"propensity"} (default) or \code{"scores"}.
#' @param ... Additional arguments passed to ggplot.
#'
#' @return A \code{ggplot} object.
#'
#' @method plot didml
#' @export
plot.didml <- function(x, type = "propensity", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Install 'ggplot2' for plotting.", call. = FALSE)

  if (type == "propensity") {
    # Need nuisance estimates for propensity plot
    if (is.null(x$nuisance)) {
      stop("Propensity plot requires DML nuisance estimates ",
           "(not available for DRDID path).", call. = FALSE)
    }
    pG_display <- if (!is.null(x$nuisance$pG_raw_original))
      x$nuisance$pG_raw_original else x$nuisance$pG_raw
    df <- data.frame(pG = pG_display)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = pG)) +
      ggplot2::geom_histogram(fill = "#012169", color = "white",
                               bins = 40, alpha = 0.8) +
      ggplot2::labs(
        x = expression(hat(e)(X) == Pr(G == 1 ~ "|" ~ X)),
        y = "Count",
        title = "Estimated propensity score distribution") +
      ggplot2::theme_minimal(base_size = 13)

    if (!is.null(x$trim_info) && x$trim_info$alpha_hat > 0) {
      p <- p + ggplot2::geom_vline(
        xintercept = 1 - x$trim_info$alpha_hat,
        linetype = "dashed", color = "#b91c1c", linewidth = 0.8
      ) +
        ggplot2::annotate(
          "text",
          x = 1 - x$trim_info$alpha_hat - 0.03,
          y = Inf, vjust = 2,
          label = paste0("alpha == ", x$trim_info$alpha_hat),
          parse = TRUE, color = "#b91c1c", size = 4)
    }
    return(p)

  } else if (type == "scores") {
    scores_list <- list()
    if (!is.null(x$wald))
      scores_list[["DML-Wald"]] <- data.frame(
        score = x$wald$scores, estimand = "DML-Wald")
    if (!is.null(x$tc))
      scores_list[["DML-TC"]] <- data.frame(
        score = x$tc$scores, estimand = "DML-TC")
    if (!is.null(x$chang))
      scores_list[["DML-Chang"]] <- data.frame(
        score = x$chang$scores, estimand = "DML-Chang")

    if (length(scores_list) == 0)
      stop("No score vectors available for plotting.", call. = FALSE)

    df <- do.call(rbind, scores_list)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = score)) +
      ggplot2::geom_histogram(fill = "#012169", color = "white",
                               bins = 50, alpha = 0.8) +
      ggplot2::facet_wrap(~estimand, scales = "free") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "#b91c1c") +
      ggplot2::labs(
        x = expression(psi[i](hat(theta))),
        y = "Count",
        title = "Score distribution at estimated theta") +
      ggplot2::theme_minimal(base_size = 13)
    return(p)
  }

  stop("Unknown plot type: '", type,
       "'. Use 'propensity' or 'scores'.", call. = FALSE)
}
