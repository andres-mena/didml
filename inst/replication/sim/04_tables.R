# =============================================================================
# 04_tables.R — Stage 4: Compile LaTeX tables from estimates (INSTANT)
#
# Loads estimates from Stage 3, produces publication-quality LaTeX tables.
# Highlights best performer per metric per sample size.
# Replaces Mean SE with SE/SD ratio (calibration diagnostic).
#
# Input:  output/simulations/estimates/estimates_n{N}.rds
# Output: output/tables/sim_main.tex
#
# Usage: Rscript code/R/sim/04_tables.R
# =============================================================================

source("code/R/sim/00_config.R")

cat("04_tables.R — Compile LaTeX tables\n\n")

# =============================================================================
# LOAD ESTIMATES
# =============================================================================
all_results <- list()
for (n_val in N_VALUES) {
  f <- estimates_file(n_val)
  if (!file.exists(f)) {
    cat(sprintf("  SKIP n=%d: %s not found\n", n_val, f))
    next
  }
  saved <- readRDS(f)
  all_results[[as.character(n_val)]] <- saved
  trim_mode <- saved$trim_used$mode
  if (!is.null(trim_mode) && trim_mode == "data-driven") {
    ok_a <- is.finite(saved$alpha_chosen)
    cat(sprintf("  Loaded n=%d (M=%d, trim=data-driven, mean_alpha=%.3f)\n",
                n_val, saved$M, mean(saved$alpha_chosen[ok_a])))
  } else {
    cat(sprintf("  Loaded n=%d (M=%d)\n", n_val, saved$M))
  }
}

if (length(all_results) == 0) stop("No estimate files found. Run 03_estimators.R first.")

# =============================================================================
# COMPUTE SUMMARY STATISTICS
# =============================================================================
summary_rows <- list()
for (n_str in names(all_results)) {
  n_val <- as.integer(n_str)
  est <- all_results[[n_str]]$estimates

  for (nm in TABLE1_NAMES) {
    e <- est[, paste0(nm, "_est")]
    s <- est[, paste0(nm, "_se")]
    cov <- est[, paste0(nm, "_cov")]
    ok <- is.finite(e) & is.finite(s) & s > 0

    mc_sd <- sd(e[ok])
    mean_se <- mean(s[ok])

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      n = n_val,
      estimator = nm,
      abs_bias = abs(mean(e[ok]) - TAU),
      sd = mc_sd,
      rmse = sqrt(mean((e[ok] - TAU)^2)),
      se_sd_ratio = mean_se / mc_sd,
      coverage = mean(cov[ok], na.rm = TRUE),
      M_valid = sum(ok),
      stringsAsFactors = FALSE
    )
  }
}
df <- do.call(rbind, summary_rows)

# =============================================================================
# HELPER: format value, bold if best
# =============================================================================
fmt_best <- function(vals, fmt, best_idx) {
  out <- sprintf(fmt, vals)
  out[best_idx] <- sprintf("\\textbf{%s}", out[best_idx])
  out
}

# =============================================================================
# BUILD LATEX TABLE
# =============================================================================
cat("\nGenerating main table with highlights...\n")

n_est <- length(TABLE1_NAMES)
tex <- c(
  sprintf("\\begin{tabular}{l%s}", paste(rep("r", n_est), collapse = "")),
  "\\toprule",
  sprintf(" & %s \\\\", paste(TABLE1_LABELS, collapse = " & ")),
  "\\midrule"
)

for (n_val in sort(unique(df$n))) {
  sub <- df[df$n == n_val, ]
  n_label <- sprintf("$n = %s$", format(n_val, big.mark = ","))

  # Best = lowest |Bias|
  best_bias <- which.min(sub$abs_bias)
  bias_str <- fmt_best(sub$abs_bias, "%.3f", best_bias)

  # Best = lowest SD
  best_sd <- which.min(sub$sd)
  sd_str <- fmt_best(sub$sd, "%.3f", best_sd)

  # Best = lowest RMSE
  best_rmse <- which.min(sub$rmse)
  rmse_str <- fmt_best(sub$rmse, "%.3f", best_rmse)

  # Best = closest to 0.95 coverage
  best_cov <- which.min(abs(sub$coverage - 0.95))
  cov_str <- fmt_best(sub$coverage, "%.2f", best_cov)

  # Best = SE/SD ratio closest to 1
  best_ratio <- which.min(abs(sub$se_sd_ratio - 1))
  ratio_str <- fmt_best(sub$se_sd_ratio, "%.2f", best_ratio)

  tex <- c(tex, sprintf("\\multicolumn{%d}{l}{%s} \\\\", n_est + 1, n_label))
  tex <- c(tex, sprintf("\\quad $|$Bias$|$ & %s \\\\", paste(bias_str, collapse = " & ")))
  tex <- c(tex, sprintf("\\quad SD & %s \\\\", paste(sd_str, collapse = " & ")))
  tex <- c(tex, sprintf("\\quad RMSE & %s \\\\", paste(rmse_str, collapse = " & ")))
  tex <- c(tex, sprintf("\\quad Coverage & %s \\\\", paste(cov_str, collapse = " & ")))
  tex <- c(tex, sprintf("\\quad $\\overline{SE}$/SD & %s \\\\", paste(ratio_str, collapse = " & ")))

  if (n_val != max(df$n)) tex <- c(tex, "\\midrule")
}

tex <- c(tex, "\\bottomrule", "\\end{tabular}")

outfile <- file.path(DIR_TABLES, "sim_main.tex")
writeLines(tex, outfile)
cat(sprintf("  Saved %s\n", outfile))

# =============================================================================
# PRINT SUMMARY
# =============================================================================
cat("\n========== SUMMARY ==========\n")
for (n_val in sort(unique(df$n))) {
  cat(sprintf("\nn = %d:\n", n_val))
  sub <- df[df$n == n_val, ]
  for (i in 1:nrow(sub)) {
    cat(sprintf("  %15s: |bias|=%6.3f  sd=%6.3f  rmse=%6.3f  cov=%.2f  SE/SD=%.2f  (%d)\n",
                sub$estimator[i], sub$abs_bias[i], sub$sd[i], sub$rmse[i],
                sub$coverage[i], sub$se_sd_ratio[i], sub$M_valid[i]))
  }
}

cat("\nDone.\n")
