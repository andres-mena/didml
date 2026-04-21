# =============================================================================
# 06_trimming_analysis.R — Trimming analysis figures for paper + appendix
#
# Produces figures that validate the optimal trimming theorem (Theorem 2)
# and demonstrate the bias-variance tradeoff.
#
# Input:  data/simulations/nuisance/nuisance_n{N}.rds (for recomputing at each alpha)
#         output/simulations/estimates/estimates_n{N}.rds (for alpha_chosen)
# Output: output/figures/sim_trimming_wald.pdf      — Bias²+Var+MSE vs α (DR-Wald)
#         output/figures/sim_trimming_tc.pdf         — Same for DR-TC (different scale)
#         output/figures/sim_trimming_alpha_dist.pdf  — Distribution of data-driven α̂
#         output/figures/sim_cost_function.pdf        — R̂(α) curve with optimal
#         output/figures/sim_kappa_validation.pdf     — κ(x) vs e(x)
#
# Usage: Rscript code/R/sim/06_trimming_analysis.R
# =============================================================================

source("code/R/sim/00_config.R")
source("code/R/04_estimators.R")
library(ggplot2)

primary_blue  <- "#012169"
primary_gold  <- "#f2a900"
accent_gray   <- "#525252"
positive_green <- "#15803d"
negative_red  <- "#b91c1c"

theme_paper <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = primary_blue, size = 13),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    )
}

ALPHA_SWEEP <- seq(0, 0.45, by = 0.01)

cat("06_trimming_analysis.R\n\n")

# =============================================================================
# FIGURE 1 & 2: BIAS² + VARIANCE + MSE vs α (separate for Wald and TC)
# =============================================================================
cat("Computing MSE decomposition across trim grid...\n")

compute_dr_at_alpha <- function(nuis, alpha) {
  Y <- nuis$Y; D <- nuis$D; G <- nuis$G; Ti <- nuis$Ti
  n <- length(Y)
  GT <- G * Ti; p_11 <- mean(GT); pT <- mean(Ti)
  if (p_11 < 0.01) return(c(NA, NA, NA, NA))

  cf <- nuis$cf
  pG <- pmin(cf$pG_raw, 1 - alpha)
  pG <- pmax(pG, 0.001)
  cf$pi_11 <- pG * pT; cf$pi_10 <- pG * (1 - pT)
  cf$pi_01 <- (1 - pG) * pT; cf$pi_00 <- (1 - pG) * (1 - pT)
  W_df <- data.frame(Y = Y, D = D, G = G, Ti = Ti)

  drw <- tryCatch(compute_dd_wald(W_df, cf), error = function(e) NULL)
  drt <- tryCatch(compute_dd_tc(W_df, cf), error = function(e) NULL)

  c(if (!is.null(drw)) drw$estimate else NA,
    if (!is.null(drw)) drw$se else NA,
    if (!is.null(drt)) drt$estimate else NA,
    if (!is.null(drt)) drt$se else NA)
}

mse_data <- list()
for (n_val in N_VALUES) {
  f <- nuisance_file(n_val)
  if (!file.exists(f)) next

  cat(sprintf("  n=%d: loading...\n", n_val))
  saved <- readRDS(f)
  all_nuis <- saved$nuisance
  M_actual <- length(all_nuis)

  for (alpha in ALPHA_SWEEP) {
    dw <- dt <- rep(NA_real_, M_actual)
    for (r in 1:M_actual) {
      nuis <- all_nuis[[r]]
      if (inherits(nuis, "error") || !is.list(nuis) || is.null(nuis$Y)) next
      res <- tryCatch(compute_dr_at_alpha(nuis, alpha), error = function(e) rep(NA, 4))
      dw[r] <- res[1]; dt[r] <- res[3]
    }

    ok_w <- is.finite(dw)
    ok_t <- is.finite(dt)
    if (sum(ok_w) > 10) {
      bias_w <- mean(dw[ok_w]) - TAU
      var_w  <- var(dw[ok_w])
      mse_data[[length(mse_data) + 1]] <- data.frame(
        n = n_val, alpha = alpha, estimator = "DR-Wald",
        bias_sq = bias_w^2, variance = var_w, mse = bias_w^2 + var_w
      )
    }
    if (sum(ok_t) > 10) {
      bias_t <- mean(dt[ok_t]) - TAU
      var_t  <- var(dt[ok_t])
      mse_data[[length(mse_data) + 1]] <- data.frame(
        n = n_val, alpha = alpha, estimator = "DR-TC",
        bias_sq = bias_t^2, variance = var_t, mse = bias_t^2 + var_t
      )
    }
  }
  cat(sprintf("  n=%d: done\n", n_val))
}

df_mse <- do.call(rbind, mse_data)

# Reshape for plotting
df_long <- rbind(
  transform(df_mse[, c("n", "alpha", "estimator", "bias_sq")],
            component = "Bias²", value = bias_sq)[, c("n", "alpha", "estimator", "component", "value")],
  transform(df_mse[, c("n", "alpha", "estimator", "variance")],
            component = "Variance", value = variance)[, c("n", "alpha", "estimator", "component", "value")],
  transform(df_mse[, c("n", "alpha", "estimator", "mse")],
            component = "MSE", value = mse)[, c("n", "alpha", "estimator", "component", "value")]
)

# Use n >= 500 to avoid extreme variance at n=100
df_long_w <- df_long[df_long$estimator == "DR-Wald" & df_long$n >= 500, ]
df_long_t <- df_long[df_long$estimator == "DR-TC" & df_long$n >= 500, ]

make_n_label <- function(n) factor(sprintf("n = %s", format(n, big.mark = ",")),
                                    levels = sprintf("n = %s", format(sort(unique(n)), big.mark = ",")))
df_long_w$n_label <- make_n_label(df_long_w$n)
df_long_t$n_label <- make_n_label(df_long_t$n)

# Mark optimal
opt_w <- do.call(rbind, lapply(split(df_mse[df_mse$estimator == "DR-Wald" & df_mse$n >= 500, ],
                                      df_mse$n[df_mse$estimator == "DR-Wald" & df_mse$n >= 500]),
                                function(x) x[which.min(x$mse), ]))
opt_t <- do.call(rbind, lapply(split(df_mse[df_mse$estimator == "DR-TC" & df_mse$n >= 500, ],
                                      df_mse$n[df_mse$estimator == "DR-TC" & df_mse$n >= 500]),
                                function(x) x[which.min(x$mse), ]))
opt_w$n_label <- make_n_label(opt_w$n)
opt_t$n_label <- make_n_label(opt_t$n)

comp_colors <- c("Bias²" = negative_red, "Variance" = primary_blue, "MSE" = accent_gray)

p_wald <- ggplot(df_long_w, aes(x = alpha, y = value, color = component)) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = opt_w, aes(xintercept = alpha), linetype = "dashed", color = accent_gray, linewidth = 0.4) +
  facet_wrap(~n_label, scales = "free_y", nrow = 1) +
  scale_color_manual(values = comp_colors) +
  labs(x = expression(alpha), y = NULL, color = NULL,
       title = "DR-Wald: Bias-variance decomposition") +
  theme_paper()

ggsave(file.path(DIR_FIGURES, "sim_trimming_wald.pdf"), p_wald,
       width = 14, height = 4, bg = "transparent")
cat("  → sim_trimming_wald.pdf\n")

p_tc <- ggplot(df_long_t, aes(x = alpha, y = value, color = component)) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = opt_t, aes(xintercept = alpha), linetype = "dashed", color = accent_gray, linewidth = 0.4) +
  facet_wrap(~n_label, scales = "free_y", nrow = 1) +
  scale_color_manual(values = comp_colors) +
  labs(x = expression(alpha), y = NULL, color = NULL,
       title = "DR-TC: Bias-variance decomposition") +
  theme_paper()

ggsave(file.path(DIR_FIGURES, "sim_trimming_tc.pdf"), p_tc,
       width = 14, height = 4, bg = "transparent")
cat("  → sim_trimming_tc.pdf\n")

# =============================================================================
# FIGURE 3: DATA-DRIVEN α̂ DISTRIBUTION
# =============================================================================
cat("Generating alpha distribution figure...\n")

alpha_data <- list()
for (n_val in N_VALUES) {
  f <- estimates_file(n_val)
  if (!file.exists(f)) next
  est <- readRDS(f)
  if (is.null(est$alpha_chosen)) next
  ok <- is.finite(est$alpha_chosen)
  alpha_data[[length(alpha_data) + 1]] <- data.frame(
    alpha = est$alpha_chosen[ok],
    n = n_val,
    n_label = sprintf("n = %s", format(n_val, big.mark = ","))
  )
}
df_alpha <- do.call(rbind, alpha_data)
df_alpha <- df_alpha[!is.na(df_alpha$n) & df_alpha$n %in% N_VALUES, ]
# Recreate n_label from n to avoid format() padding mismatch between scalar and vector
df_alpha$n_label <- sprintf("n = %s", formatC(df_alpha$n, format = "d", big.mark = ","))
df_alpha$n_label <- factor(df_alpha$n_label,
                           levels = sprintf("n = %s", formatC(sort(unique(df_alpha$n)), format = "d", big.mark = ",")))

p_alpha <- ggplot(df_alpha, aes(x = alpha)) +
  geom_histogram(fill = primary_blue, color = "white", bins = 30, alpha = 0.8) +
  facet_wrap(~n_label, scales = "free_y", nrow = 1) +
  labs(x = expression(hat(alpha)), y = "Count",
       title = expression("Distribution of data-driven " * hat(alpha) * " across replications")) +
  theme_paper()

ggsave(file.path(DIR_FIGURES, "sim_trimming_alpha_dist.pdf"), p_alpha,
       width = 14, height = 3.5, bg = "transparent")
cat("  → sim_trimming_alpha_dist.pdf\n")

# =============================================================================
# FIGURE 4: R̂(α) CURVE (representative replication at n=2000)
# =============================================================================
cat("Generating R_hat(alpha) curve...\n")

n_demo <- 2000
f_demo <- nuisance_file(n_demo)
if (file.exists(f_demo)) {
  saved_demo <- readRDS(f_demo)
  nuis_demo <- saved_demo$nuisance[[1]]  # rep 1

  pG_raw <- nuis_demo$cf$pG_raw
  DID_D_hat <- nuis_demo$D - nuis_demo$cf$m_D_10 - nuis_demo$cf$m_D_01 + nuis_demo$cf$m_D_00

  # Estimate g_hat
  g_fit <- loess(DID_D_hat ~ pG_raw, span = 0.3, degree = 1)
  g_hat <- pmax(predict(g_fit), 0.001)
  g_hat[is.na(g_hat)] <- mean(DID_D_hat)

  R_curve <- data.frame(alpha = ALPHA_SWEEP, R_hat = NA_real_)
  for (i in seq_along(ALPHA_SWEEP)) {
    a <- ALPHA_SWEEP[i]
    keep <- pG_raw <= (1 - a)
    if (sum(keep) < 10) next
    num <- sum(pG_raw[keep] / (1 - pG_raw[keep]))
    den <- sum(pG_raw[keep] * g_hat[keep])^2
    if (den > 1e-10) R_curve$R_hat[i] <- num / den
  }

  opt_idx <- which.min(R_curve$R_hat)

  p_rcurve <- ggplot(R_curve[is.finite(R_curve$R_hat), ], aes(x = alpha, y = R_hat)) +
    geom_line(color = primary_blue, linewidth = 1) +
    geom_vline(xintercept = R_curve$alpha[opt_idx], linetype = "dashed", color = negative_red) +
    annotate("text", x = R_curve$alpha[opt_idx] + 0.02,
             y = max(R_curve$R_hat[is.finite(R_curve$R_hat)]) * 0.9,
             label = sprintf("hat(alpha) == %.2f", R_curve$alpha[opt_idx]),
             parse = TRUE, color = negative_red, size = 4) +
    labs(x = expression(alpha), y = expression(hat(R)(alpha)),
         title = expression("Variance objective " * hat(R)(alpha) * " (Remark 2, n = 2,000, rep 1)")) +
    theme_paper()

  ggsave(file.path(DIR_FIGURES, "sim_cost_function.pdf"), p_rcurve,
         width = 7, height = 4, bg = "transparent")
  cat("  → sim_cost_function.pdf\n")
}

# =============================================================================
# FIGURE 5: κ(x) vs e(x) VALIDATION
# =============================================================================
cat("Generating kappa validation figure...\n")

if (file.exists(f_demo)) {
  # Compute empirical kappa for a representative replication
  # κ^wald(x) ∝ e/(1-e) / (e * DID_D) = 1/((1-e)*DID_D)
  # Under homoskedasticity: κ ∝ 1/((1-e)*g(e))
  n_obs <- length(pG_raw)

  kappa_empirical <- (pG_raw / (1 - pmax(pG_raw, 0.01))) / (pG_raw * pmax(g_hat, 0.001))

  # Theoretical curve under homoskedasticity with constant g
  e_grid <- seq(0.05, 0.95, by = 0.01)
  g_bar <- mean(DID_D_hat)
  kappa_theory_const <- 1 / ((1 - e_grid) * g_bar)
  kappa_theory_general <- 1 / ((1 - e_grid) * pmax(predict(g_fit, newdata = data.frame(pG_raw = e_grid)), 0.001))
  kappa_theory_general[is.na(kappa_theory_general)] <- NA

  df_kappa <- data.frame(
    e = pG_raw,
    kappa = kappa_empirical
  )
  # Trim extreme kappa for visualization
  q99 <- quantile(kappa_empirical, 0.99, na.rm = TRUE)
  df_kappa <- df_kappa[is.finite(df_kappa$kappa) & df_kappa$kappa < q99, ]

  df_theory <- rbind(
    data.frame(e = e_grid, kappa = kappa_theory_const, type = "Cor 2.2 (constant g)"),
    data.frame(e = e_grid, kappa = kappa_theory_general, type = "Cor 2.1 (estimated g)")
  )

  p_kappa <- ggplot() +
    geom_point(data = df_kappa, aes(x = e, y = kappa), alpha = 0.15, size = 0.8, color = accent_gray) +
    geom_line(data = df_theory, aes(x = e, y = kappa, color = type, linetype = type), linewidth = 1) +
    scale_color_manual(values = c("Cor 2.2 (constant g)" = primary_gold,
                                   "Cor 2.1 (estimated g)" = primary_blue)) +
    scale_linetype_manual(values = c("Cor 2.2 (constant g)" = "dashed",
                                      "Cor 2.1 (estimated g)" = "solid")) +
    labs(x = expression(hat(e)(X)), y = expression(hat(kappa)(x)),
         color = NULL, linetype = NULL,
         title = expression("Identification-adjusted cost " * kappa(x) * " vs propensity (n = 2,000, rep 1)")) +
    theme_paper()

  ggsave(file.path(DIR_FIGURES, "sim_kappa_validation.pdf"), p_kappa,
         width = 7, height = 4.5, bg = "transparent")
  cat("  → sim_kappa_validation.pdf\n")
}

# Save MSE grid for other scripts
saveRDS(list(mse = df_mse, optimal_wald = opt_w, optimal_tc = opt_t,
             timestamp = Sys.time()),
        file.path(DIR_ESTIMATES, "trimming_mse_grid.rds"))
cat(sprintf("  → %s\n", file.path(DIR_ESTIMATES, "trimming_mse_grid.rds")))

cat("\n=== Trimming analysis complete ===\n")
