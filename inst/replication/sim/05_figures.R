# =============================================================================
# 05_figures.R — Stage 5: Generate publication figures from estimates (INSTANT)
#
# Loads estimates from Stage 3, produces figures for the paper.
# No estimation, no ML, no heavy computation.
#
# Input:  output/simulations/estimates/estimates_n{N}.rds
# Output: output/figures/sim_density.pdf       — kernel density of estimates
#         output/figures/sim_bias_variance.pdf  — bias-variance tradeoff
#         output/figures/sim_coverage.pdf       — coverage rates
#
# Usage:
#   Rscript code/R/sim/05_figures.R
#
# =============================================================================

source("code/R/sim/00_config.R")
library(ggplot2)

# Visual identity
primary_blue  <- "#012169"
primary_gold  <- "#f2a900"
accent_gray   <- "#525252"
positive_green <- "#15803d"
negative_red  <- "#b91c1c"

theme_paper <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = primary_blue),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    )
}

# Estimator colors (7 estimators)
est_colors <- c(
  "wald" = accent_gray,
  "waldx_lasso" = "#4e79a7",
  "waldx_sieve" = "#76b7b2",
  "tcx_lasso" = primary_gold,
  "tcx_sieve" = "#edc948",
  "drwald" = primary_blue,
  "drtc" = negative_red
)

est_labels_short <- c(
  "wald" = "Wald",
  "waldx_lasso" = "Wald_X(L)",
  "waldx_sieve" = "Wald_X(S)",
  "tcx_lasso" = "TC_X(L)",
  "tcx_sieve" = "TC_X(S)",
  "drwald" = "DR-Wald",
  "drtc" = "DR-TC"
)

# =============================================================================
# LOAD ESTIMATES
# =============================================================================
cat("05_figures.R — Generate simulation figures\n\n")

all_results <- list()
for (n_val in N_VALUES) {
  f <- estimates_file(n_val)
  if (!file.exists(f)) next
  all_results[[as.character(n_val)]] <- readRDS(f)
}
if (length(all_results) == 0) stop("No estimate files found. Run 03_estimators.R first.")

# =============================================================================
# FIGURE 1: DENSITY OF ESTIMATES (faceted by n)
# =============================================================================
cat("  Figure 1: density of estimates...\n")

# Focus on the 4 key estimators for clarity
key_estimators <- c("wald", "waldx_lasso", "drwald", "drtc")

density_data <- list()
for (n_str in names(all_results)) {
  est <- all_results[[n_str]]$estimates
  n_val <- as.integer(n_str)
  for (nm in key_estimators) {
    e <- est[, paste0(nm, "_est")]
    ok <- is.finite(e)
    if (sum(ok) < 10) next
    density_data[[length(density_data) + 1]] <- data.frame(
      estimate = e[ok],
      estimator = nm,
      n = n_val,
      n_label = sprintf("n = %s", format(n_val, big.mark = ","))
    )
  }
}
df_dens <- do.call(rbind, density_data)
df_dens$n_label <- factor(df_dens$n_label,
                          levels = sprintf("n = %s", format(sort(unique(df_dens$n)), big.mark = ",")))
df_dens$estimator <- factor(df_dens$estimator, levels = key_estimators)

# Trim extreme outliers for visualization
for (n_val in unique(df_dens$n)) {
  idx <- df_dens$n == n_val
  q <- quantile(df_dens$estimate[idx], c(0.01, 0.99), na.rm = TRUE)
  df_dens <- df_dens[!(idx & (df_dens$estimate < q[1] | df_dens$estimate > q[2])), ]
}

p1 <- ggplot(df_dens, aes(x = estimate, fill = estimator, color = estimator)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  geom_vline(xintercept = TAU, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~n_label, scales = "free", nrow = 1) +
  scale_fill_manual(values = est_colors[key_estimators], labels = est_labels_short[key_estimators]) +
  scale_color_manual(values = est_colors[key_estimators], labels = est_labels_short[key_estimators]) +
  labs(x = "Estimate", y = "Density", fill = NULL, color = NULL) +
  theme_paper()

ggsave(file.path(DIR_FIGURES, "sim_density.pdf"), p1,
       width = 14, height = 4, bg = "transparent")
cat(sprintf("    → %s\n", file.path(DIR_FIGURES, "sim_density.pdf")))

# =============================================================================
# FIGURE 2: BIAS vs SD SCATTERPLOT (all 7 estimators, faceted by n)
# =============================================================================
cat("  Figure 2: bias vs SD...\n")

summary_data <- list()
for (n_str in names(all_results)) {
  est <- all_results[[n_str]]$estimates
  n_val <- as.integer(n_str)
  for (nm in ESTIMATOR_NAMES) {
    e <- est[, paste0(nm, "_est")]
    ok <- is.finite(e)
    if (sum(ok) < 10) next
    summary_data[[length(summary_data) + 1]] <- data.frame(
      bias = abs(mean(e[ok]) - TAU),
      sd = sd(e[ok]),
      estimator = nm,
      n = n_val,
      n_label = sprintf("n = %s", format(n_val, big.mark = ","))
    )
  }
}
df_sum <- do.call(rbind, summary_data)
df_sum$n_label <- factor(df_sum$n_label,
                         levels = sprintf("n = %s", format(sort(unique(df_sum$n)), big.mark = ",")))

p2 <- ggplot(df_sum, aes(x = bias, y = sd, color = estimator, shape = estimator)) +
  geom_point(size = 3.5) +
  facet_wrap(~n_label, scales = "free", nrow = 1) +
  scale_color_manual(values = est_colors, labels = est_labels_short) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 16, 17), labels = est_labels_short) +
  labs(x = "|Bias|", y = "Standard Deviation", color = NULL, shape = NULL) +
  theme_paper()

ggsave(file.path(DIR_FIGURES, "sim_bias_variance.pdf"), p2,
       width = 14, height = 4, bg = "transparent")
cat(sprintf("    → %s\n", file.path(DIR_FIGURES, "sim_bias_variance.pdf")))

# =============================================================================
# FIGURE 3: COVERAGE RATES (bar chart, faceted by n)
# =============================================================================
cat("  Figure 3: coverage rates...\n")

cov_data <- list()
for (n_str in names(all_results)) {
  est <- all_results[[n_str]]$estimates
  n_val <- as.integer(n_str)
  for (nm in ESTIMATOR_NAMES) {
    cov <- est[, paste0(nm, "_cov")]
    ok <- is.finite(cov)
    if (sum(ok) < 10) next
    cov_data[[length(cov_data) + 1]] <- data.frame(
      coverage = mean(cov[ok]),
      estimator = nm,
      n = n_val,
      n_label = sprintf("n = %s", format(n_val, big.mark = ","))
    )
  }
}
df_cov <- do.call(rbind, cov_data)
df_cov$n_label <- factor(df_cov$n_label,
                         levels = sprintf("n = %s", format(sort(unique(df_cov$n)), big.mark = ",")))
df_cov$estimator <- factor(df_cov$estimator, levels = ESTIMATOR_NAMES)

p3 <- ggplot(df_cov, aes(x = estimator, y = coverage, fill = estimator)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  facet_wrap(~n_label, nrow = 1) +
  scale_fill_manual(values = est_colors) +
  scale_x_discrete(labels = est_labels_short) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Coverage (95% CI)", fill = NULL) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none")

ggsave(file.path(DIR_FIGURES, "sim_coverage.pdf"), p3,
       width = 14, height = 4, bg = "transparent")
cat(sprintf("    → %s\n", file.path(DIR_FIGURES, "sim_coverage.pdf")))

cat("\n=== Stage 5 complete ===\n")
