# =============================================================================
# 08_figures_sensitivity.R — regenerate sensitivity figures with fixed y-axis
#
# Reads J01 (primary, Lasso) and K01 (high-school, Lasso) estimates from the
# empirical_ext exploration. Plots DML-Wald and DML-TC point estimates over
# alpha in [0, 0.20] with a common y-axis [0, 2.5] for side-by-side comparison.
#
# Output:
#   output/figures/empirical_sensitivity_primary.pdf
#   output/figures/empirical_sensitivity_highschool.pdf
# =============================================================================

source("code/R/sim/00_config.R")
library(ggplot2)

primary_blue  <- "#012169"
accent_gray   <- "#525252"
negative_red  <- "#b91c1c"

theme_paper <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", color = primary_blue, size = 13),
      plot.subtitle = element_text(color = accent_gray, size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 11),
      axis.text  = element_text(size = 10)
    )
}

build_sensitivity <- function(est_path, row_path, title) {
  e   <- readRDS(est_path)
  row <- read.csv(row_path, stringsAsFactors = FALSE)

  alpha_hat <- e$alpha_hat
  dr_grid <- e$dr_grid
  df <- do.call(rbind, lapply(dr_grid, function(r) {
    data.frame(
      alpha = c(r$alpha, r$alpha),
      estimator = c("DML-Wald", "DML-TC"),
      estimate  = c(r$drwald, r$drtc),
      se        = c(r$drwald_se, r$drtc_se)
    )
  }))
  df <- df[!is.na(df$estimate), ]

  ggplot(df, aes(x = alpha, y = estimate, color = estimator)) +
    geom_ribbon(aes(ymin = estimate - 1.96 * se,
                    ymax = estimate + 1.96 * se,
                    fill = estimator),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_vline(xintercept = alpha_hat, linetype = "dashed",
               color = accent_gray, linewidth = 0.5) +
    scale_color_manual(values = c("DML-Wald" = primary_blue,
                                  "DML-TC"   = negative_red)) +
    scale_fill_manual(values = c("DML-Wald" = primary_blue,
                                 "DML-TC"   = negative_red)) +
    scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5),
                       expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 0.20),
                       breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                       expand = c(0.01, 0.01)) +
    labs(x = expression(alpha ~ "(trimming threshold)"),
         y = "Point estimate (log points per complier)",
         color = NULL, fill = NULL,
         title = title) +
    theme_paper()
}

# Primary completion (J01)
p_primary <- build_sensitivity(
  est_path = "explorations/empirical_ext/data/J01_estimates.rds",
  row_path = "explorations/empirical_ext/summary/J01_row.csv",
  title = "Primary completion"
)
ggsave(file.path(DIR_FIGURES, "empirical_sensitivity_primary.pdf"),
       p_primary, width = 7, height = 5, bg = "transparent")
message(sprintf("Saved: %s",
                file.path(DIR_FIGURES, "empirical_sensitivity_primary.pdf")))

# High-school completion (K01)
p_highschool <- build_sensitivity(
  est_path = "explorations/empirical_ext/data/K01_estimates.rds",
  row_path = "explorations/empirical_ext/summary/K01_row.csv",
  title = "High-school completion"
)
ggsave(file.path(DIR_FIGURES, "empirical_sensitivity_highschool.pdf"),
       p_highschool, width = 7, height = 5, bg = "transparent")
message(sprintf("Saved: %s",
                file.path(DIR_FIGURES, "empirical_sensitivity_highschool.pdf")))
