# =============================================================================
# trim_rule_summary.R --- Summarize data-driven one-sided vs Crump symmetric
#
# Loads saved trim_comparison_n{N}.rds (from compare_trimming.R) and the raw
# nuisance data, and produces a table comparing:
#
#   A. Data-driven one-sided   (pG <- pmin(pG_raw, 1 - hat_alpha))
#   B. Crump symmetric [0.10, 0.90]   (pG <- pmax(pmin(pG_raw, 0.90), 0.10))
#
# Metrics: Bias, SD, RMSE, Coverage for DR-Wald and DR-TC,
# plus effective sample size (share of observations with pG in (0.10, 0.90)
# versus share kept under the one-sided rule).
#
# Output: output/tables/sim_trim_rule_summary.tex  +  console print-out.
# =============================================================================

source("code/R/sim/00_config.R")

N_VALUES_AVAIL <- c(100, 500, 1000, 2000)  # n=5000 not yet in trim_comparison
results <- list()

cat("=================================================================\n")
cat("Trim-rule comparison: data-driven one-sided vs Crump symmetric\n")
cat("=================================================================\n\n")

for (n in N_VALUES_AVAIL) {
  comp_file <- file.path(DIR_ESTIMATES, sprintf("trim_comparison_n%d.rds", n))
  nuis_file <- nuisance_file(n)

  if (!file.exists(comp_file) || !file.exists(nuis_file)) {
    cat(sprintf("SKIP n=%d (missing files)\n", n))
    next
  }

  comp <- readRDS(comp_file)
  nuis_all <- readRDS(nuis_file)$nuisance
  M <- comp$M

  # ---- metrics from saved estimates ----
  metric_for <- function(est_mat, se_mat, col) {
    e <- est_mat[, col]; s <- se_mat[, col]
    ok <- is.finite(e) & is.finite(s) & s > 0
    if (sum(ok) < 10) return(c(NA, NA, NA, NA))
    bias <- mean(e[ok]) - comp$tau
    sdv  <- sd(e[ok])
    rmse <- sqrt(mean((e[ok] - comp$tau)^2))
    cov  <- mean((e[ok] - 1.96 * s[ok] <= comp$tau) &
                 (comp$tau <= e[ok] + 1.96 * s[ok]))
    c(bias, sdv, rmse, cov)
  }

  w_dd  <- metric_for(comp$drwald_est, comp$drwald_se, "data_1s")
  w_sym <- metric_for(comp$drwald_est, comp$drwald_se, "fixed_010")
  t_dd  <- metric_for(comp$drtc_est,   comp$drtc_se,   "data_1s")
  t_sym <- metric_for(comp$drtc_est,   comp$drtc_se,   "fixed_010")

  # ---- effective sample share per replication ----
  share_kept_dd  <- numeric(M)  # fraction with pG_raw <= 1 - hat_alpha
  share_free_sym <- numeric(M)  # fraction with pG_raw in (0.10, 0.90)
  hat_alphas     <- numeric(M)

  for (r in seq_len(M)) {
    nuis <- nuis_all[[r]]
    if (inherits(nuis, "error") || is.null(nuis$cf$pG_raw)) {
      share_kept_dd[r] <- share_free_sym[r] <- hat_alphas[r] <- NA
      next
    }
    p <- nuis$cf$pG_raw
    a <- comp$alpha_chosen[r, "data_1s"]
    hat_alphas[r]     <- a
    share_kept_dd[r]  <- mean(p <= (1 - a))
    share_free_sym[r] <- mean(p > 0.10 & p < 0.90)
  }

  mean_kept_dd  <- mean(share_kept_dd,  na.rm = TRUE)
  mean_free_sym <- mean(share_free_sym, na.rm = TRUE)
  mean_alpha    <- mean(hat_alphas,     na.rm = TRUE)

  results[[as.character(n)]] <- list(
    n = n, w_dd = w_dd, w_sym = w_sym, t_dd = t_dd, t_sym = t_sym,
    mean_kept_dd = mean_kept_dd, mean_free_sym = mean_free_sym,
    mean_alpha = mean_alpha
  )

  # ---- console print ----
  cat(sprintf("---- n = %d (M = %d) ----\n", n, M))
  cat(sprintf("  Data-driven hat_alpha: mean = %.3f\n", mean_alpha))
  cat(sprintf("  Share of obs with uncapped weight:\n"))
  cat(sprintf("     Data-driven one-sided: %.1f%% (pG_raw <= 1 - hat_alpha)\n",
              100 * mean_kept_dd))
  cat(sprintf("     Crump symmetric [.1, .9]: %.1f%% (pG_raw in (.1, .9))\n\n",
              100 * mean_free_sym))

  cat(sprintf("  %-28s %8s %8s %8s %8s\n", "Rule", "Bias", "SD", "RMSE", "Cov"))
  cat(sprintf("  %-28s %8.3f %8.3f %8.3f %8.2f\n",
              "DR-Wald | Data-driven 1-sided", w_dd[1], w_dd[2], w_dd[3], w_dd[4]))
  cat(sprintf("  %-28s %8.3f %8.3f %8.3f %8.2f\n",
              "DR-Wald | Symmetric [.1, .9]",   w_sym[1], w_sym[2], w_sym[3], w_sym[4]))
  cat(sprintf("  %-28s %8.3f %8.3f %8.3f %8.2f\n",
              "DR-TC   | Data-driven 1-sided", t_dd[1], t_dd[2], t_dd[3], t_dd[4]))
  cat(sprintf("  %-28s %8.3f %8.3f %8.3f %8.2f\n",
              "DR-TC   | Symmetric [.1, .9]",   t_sym[1], t_sym[2], t_sym[3], t_sym[4]))
  cat("\n")
}

# =============================================================================
# LaTeX TABLE
# =============================================================================
tex <- character(0)
tex <- c(tex, "\\begin{tabular}{cl rrrr rr}")
tex <- c(tex, "\\toprule")
tex <- c(tex, " & & \\multicolumn{4}{c}{Estimator} & \\multicolumn{2}{c}{Share uncapped} \\\\")
tex <- c(tex, "\\cmidrule(lr){3-6} \\cmidrule(lr){7-8}")
tex <- c(tex, "$n$ & Trimming rule & Bias & SD & RMSE & Cov. & One-sided & Symmetric \\\\")
tex <- c(tex, "\\midrule")

for (key in names(results)) {
  r <- results[[key]]
  # DR-Wald rows
  tex <- c(tex, sprintf(
    "\\multirow{4}{*}{%d} & DR-Wald, data-driven 1-sided  & %.3f & %.3f & %.3f & %.2f & \\multirow{4}{*}{%.1f\\%%} & \\multirow{4}{*}{%.1f\\%%} \\\\",
    r$n, r$w_dd[1], r$w_dd[2], r$w_dd[3], r$w_dd[4],
    100 * r$mean_kept_dd, 100 * r$mean_free_sym))
  tex <- c(tex, sprintf(
    " & DR-Wald, symmetric $[0.10, 0.90]$ & %.3f & %.3f & %.3f & %.2f & & \\\\",
    r$w_sym[1], r$w_sym[2], r$w_sym[3], r$w_sym[4]))
  # DR-TC rows
  tex <- c(tex, sprintf(
    " & DR-TC, data-driven 1-sided    & %.3f & %.3f & %.3f & %.2f & & \\\\",
    r$t_dd[1], r$t_dd[2], r$t_dd[3], r$t_dd[4]))
  tex <- c(tex, sprintf(
    " & DR-TC, symmetric $[0.10, 0.90]$   & %.3f & %.3f & %.3f & %.2f & & \\\\",
    r$t_sym[1], r$t_sym[2], r$t_sym[3], r$t_sym[4]))
  tex <- c(tex, "\\midrule")
}

tex[length(tex)] <- "\\bottomrule"
tex <- c(tex, "\\end{tabular}")

out_tex <- file.path("output/tables", "sim_trim_rule_summary.tex")
writeLines(tex, out_tex)
cat(sprintf("LaTeX table written: %s\n", out_tex))

# Also save RDS summary
saveRDS(results, file.path("output/simulations", "trim_rule_summary.rds"))
cat(sprintf("RDS summary written: %s\n",
            file.path("output/simulations", "trim_rule_summary.rds")))
