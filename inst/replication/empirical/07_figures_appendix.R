# =============================================================================
# 07_figures_appendix.R — regenerate appendix propensity figures with p=45
#
# Reads Lasso (J01) and Ridge (J02) cross-fitted nuisances from the
# empirical_ext exploration (p=45 individual-level dictionary, N=22,414).
# Plots Pr(G=1|X) by (G,T) cell for each learner.
#
# Output:
#   output/figures/empirical_propensity_lasso.pdf
#   output/figures/empirical_propensity_ridge.pdf
# =============================================================================

source("code/R/sim/00_config.R")
library(ggplot2)

primary_blue  <- "#012169"
primary_gold  <- "#f2a900"
accent_gray   <- "#525252"
positive_green <- "#15803d"
negative_red  <- "#b91c1c"

theme_paper <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", color = primary_blue, size = 14),
      plot.subtitle = element_text(color = accent_gray, size = 11),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
}

build_plot <- function(pG, G, Ti, title, alpha_cutoff) {
  cell <- ifelse(G == 1 & Ti == 1, "G=1, T=1 (treated post)",
          ifelse(G == 1 & Ti == 0, "G=1, T=0 (treated pre)",
          ifelse(G == 0 & Ti == 1, "G=0, T=1 (control post)",
                                    "G=0, T=0 (control pre)")))
  cell <- factor(cell, levels = c("G=0, T=0 (control pre)",
                                  "G=0, T=1 (control post)",
                                  "G=1, T=0 (treated pre)",
                                  "G=1, T=1 (treated post)"))
  df <- data.frame(pG = pG, cell = cell)

  ggplot(df, aes(x = pG)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40,
                   fill = primary_blue, alpha = 0.7,
                   color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 1 - alpha_cutoff, linetype = "dashed",
               color = negative_red, linewidth = 0.6) +
    facet_wrap(~ cell, ncol = 2, scales = "free_y") +
    labs(title = title,
         subtitle = sprintf("N = %s, p = 45 individual-level dictionary, $\\hat\\alpha$-rule cutoff = %.2f",
                            format(length(pG), big.mark = ","), alpha_cutoff),
         x = expression(hat(e)(X) == widehat(Pr)(G == 1 ~ "|" ~ X)),
         y = "Density") +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    theme_paper()
}

# Lasso (J01)
j1_nuis <- readRDS("explorations/empirical_ext/data/J01_nuisance.rds")
j1_dat  <- readRDS("explorations/empirical_ext/data/J01_analysis.rds")
j1_row  <- read.csv("explorations/empirical_ext/summary/J01_row.csv",
                    stringsAsFactors = FALSE)

alpha_lasso <- j1_row$alpha_hat[1]
p_lasso <- build_plot(j1_nuis$cf$pG_raw, j1_dat$G, j1_dat$Ti,
                      title = "Lasso logit first step", alpha_cutoff = alpha_lasso)

ggsave(file.path(DIR_FIGURES, "empirical_propensity_lasso.pdf"), p_lasso,
       width = 8, height = 6, bg = "transparent")
message(sprintf("Saved: %s (alpha_hat = %.3f)",
                file.path(DIR_FIGURES, "empirical_propensity_lasso.pdf"),
                alpha_lasso))

# Ridge (J02)
j2_nuis <- readRDS("explorations/empirical_ext/data/J02_nuisance.rds")
j2_dat  <- readRDS("explorations/empirical_ext/data/J02_analysis.rds")
j2_row  <- read.csv("explorations/empirical_ext/summary/J02_row.csv",
                    stringsAsFactors = FALSE)

alpha_ridge <- j2_row$alpha_hat[1]
p_ridge <- build_plot(j2_nuis$cf$pG_raw, j2_dat$G, j2_dat$Ti,
                      title = "Ridge logit first step", alpha_cutoff = alpha_ridge)

ggsave(file.path(DIR_FIGURES, "empirical_propensity_ridge.pdf"), p_ridge,
       width = 8, height = 6, bg = "transparent")
message(sprintf("Saved: %s (alpha_hat = %.3f)",
                file.path(DIR_FIGURES, "empirical_propensity_ridge.pdf"),
                alpha_ridge))

# Diagnostic summary
cat("\n=== Diagnostics ===\n")
cat(sprintf("Lasso pG: min=%.3f med=%.3f mean=%.3f max=%.3f  sd=%.3f\n",
            min(j1_nuis$cf$pG_raw), median(j1_nuis$cf$pG_raw),
            mean(j1_nuis$cf$pG_raw), max(j1_nuis$cf$pG_raw),
            sd(j1_nuis$cf$pG_raw)))
cat(sprintf("Ridge pG: min=%.3f med=%.3f mean=%.3f max=%.3f  sd=%.3f\n",
            min(j2_nuis$cf$pG_raw), median(j2_nuis$cf$pG_raw),
            mean(j2_nuis$cf$pG_raw), max(j2_nuis$cf$pG_raw),
            sd(j2_nuis$cf$pG_raw)))
cat(sprintf("Lasso share below 0.1: %.3f\n",
            mean(j1_nuis$cf$pG_raw < 0.1)))
cat(sprintf("Lasso share above 0.9: %.3f\n",
            mean(j1_nuis$cf$pG_raw > 0.9)))
cat(sprintf("Ridge share below 0.1: %.3f\n",
            mean(j2_nuis$cf$pG_raw < 0.1)))
cat(sprintf("Ridge share above 0.9: %.3f\n",
            mean(j2_nuis$cf$pG_raw > 0.9)))
