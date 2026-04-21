# =============================================================================
# 07_ml_comparison_table.R — ML method comparison table for appendix
#
# Compares DR-Wald and DR-TC across Lasso, RF, NN first-step methods.
# Uses estimates from Stage 3 (must run 03_estimators.R --method {rf,nn} first).
#
# Output: output/tables/sim_ml_comparison.tex
# =============================================================================

source("code/R/sim/00_config.R")

methods <- c("lasso", "rf", "nn")

# Use all sample sizes that have estimates for ALL methods
n_vals <- N_VALUES[sapply(N_VALUES, function(n) {
  all(sapply(methods, function(m) file.exists(estimates_file(n, m))))
})]
cat(sprintf("  Sample sizes with all methods: %s\n", paste(n_vals, collapse = ", ")))

cat("=== ML Method Comparison ===\n\n")

results <- list()
for (method in methods) {
  for (n_val in n_vals) {
    f <- estimates_file(n_val, method)
    if (!file.exists(f)) {
      cat(sprintf("  SKIP: %s\n", f))
      next
    }
    est <- readRDS(f)
    e <- est$estimates
    for (nm in c("drwald", "drtc")) {
      ev <- e[, paste0(nm, "_est")]
      sv <- e[, paste0(nm, "_se")]
      cv <- e[, paste0(nm, "_cov")]
      ok <- is.finite(ev) & is.finite(sv) & sv > 0
      results[[length(results) + 1]] <- data.frame(
        method = method, n = n_val, estimator = nm,
        bias = round(mean(ev[ok]) - TAU, 3),
        sd = round(sd(ev[ok]), 3),
        rmse = round(sqrt(mean((ev[ok] - TAU)^2)), 3),
        coverage = round(mean(cv[ok], na.rm = TRUE), 2),
        valid = sum(ok))
    }
  }
}

df <- do.call(rbind, results)

# LaTeX table
method_labels <- c(lasso = "Lasso", rf = "Random Forest", nn = "Neural Net")

tex <- c(
  "\\begin{tabular}{llrrrr}",
  "\\toprule",
  " & & \\multicolumn{2}{c}{DR-Wald} & \\multicolumn{2}{c}{DR-TC} \\\\",
  "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}",
  "$N$ & Method & RMSE & Cov. & RMSE & Cov. \\\\",
  "\\midrule"
)

for (nv in n_vals) {
  for (i in seq_along(methods)) {
    m <- methods[i]
    dw <- df[df$n == nv & df$method == m & df$estimator == "drwald", ]
    dt <- df[df$n == nv & df$method == m & df$estimator == "drtc", ]
    if (nrow(dw) == 0 || nrow(dt) == 0) next
    n_label <- if (i == 1) format(nv, big.mark = ",") else ""
    tex <- c(tex, sprintf("%s & %s & %.3f & %.2f & %.3f & %.2f \\\\",
                          n_label, method_labels[m], dw$rmse, dw$coverage,
                          dt$rmse, dt$coverage))
  }
  if (nv != max(n_vals)) tex <- c(tex, "\\midrule")
}

tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(DIR_TABLES, "sim_ml_comparison.tex"))
cat(sprintf("Saved: %s\n\n", file.path(DIR_TABLES, "sim_ml_comparison.tex")))

cat("Summary:\n")
print(df[, c("n", "method", "estimator", "bias", "sd", "rmse", "coverage")])
