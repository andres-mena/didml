# =============================================================================
# 04_tables.R -- Stage 4: Generate publication LaTeX tables (INSTANT)
#
# Produces three tables from the Stage 3 estimates:
#   1. empirical_main_v2.tex    -- 9-column main results (with clustered SEs)
#   2. empirical_comparison.tex -- analysis of estimator differences
#   3. empirical_placebo.tex    -- pre-trend / placebo test results
#
# Input:  output/simulations/empirical_estimates.rds
# Output: output/tables/empirical_main_v2.tex
#         output/tables/empirical_comparison.tex
#         output/tables/empirical_placebo.tex
#
# Author:  Andres Mena (Brown University)
# Project: fuzzyDDML
# =============================================================================

source("code/R/sim/00_config.R")

# =============================================================================
# 1. LOAD ESTIMATES
# =============================================================================

message("=== 04_tables.R: Generate empirical tables ===\n")

args <- commandArgs(trailingOnly = TRUE)
use_duflo <- "--duflo" %in% args

infile <- sprintf("output/simulations/empirical_estimates_%s.rds",
                  if (use_duflo) "duflo" else "medicaid")
if (!file.exists(infile)) stop("Run 03_estimators.R first: ", infile)
est <- readRDS(infile)
message(sprintf("  Loaded: %s", infile))

# Extract results
wald   <- est$wald
waldx  <- est$waldx
tcx    <- est$tcx
drw_o  <- est$drwald$optimal
drw_n  <- est$drwald$notrim
drw_f  <- est$drwald$fixed
drt_o  <- est$drtc$optimal
drt_n  <- est$drtc$notrim
drt_f  <- est$drtc$fixed

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

fmt_est <- function(x) sprintf("%.3f", x)
fmt_se  <- function(x) sprintf("(%.3f)", x)
fmt_ci  <- function(lo, hi) sprintf("[%.2f, %.2f]", lo, hi)

# Stars for individual estimates (two-sided test H0: theta = 0)
add_stars <- function(est_val, se_val) {
  if (is.na(est_val) || is.na(se_val) || se_val <= 0) return(fmt_est(est_val))
  z <- abs(est_val / se_val)
  s <- if (z > 2.576) "***" else if (z > 1.960) "**" else if (z > 1.645) "*" else ""
  sprintf("%.3f%s", est_val, s)
}

# =============================================================================
# 3. TABLE 1: MAIN RESULTS (9 columns)
# =============================================================================

message("Generating empirical_main_v2.tex...")

# Column order: Wald, Wald_X, DR-Wald(symmetric), DR-Wald(optimal),
#               TC_X, DR-TC(symmetric), DR-TC(optimal)
results <- list(wald, waldx, drw_f, drw_o, tcx, drt_f, drt_o)

# Build estimate row with stars
est_row <- sapply(results, function(r) add_stars(r$estimate, r$se))

# Build SE row
se_row <- sapply(results, function(r) fmt_se(r$se))

# Build clustered SE row (state-level)
se_clust_row <- sapply(results, function(r) {
  if (!is.null(r$se_clustered) && !is.na(r$se_clustered)) {
    fmt_se(r$se_clustered)
  } else {
    "---"
  }
})

# Build CI row
ci_row <- sapply(results, function(r) fmt_ci(r$ci_lower, r$ci_upper))

# Metadata rows
did_Y_str <- fmt_est(est$did_Y)
did_D_str <- fmt_est(est$did_D)
N_str     <- format(est$N, big.mark = ",")
p_str     <- format(est$p, big.mark = ",")
states_str <- if (!is.null(est$n_states)) as.character(est$n_states) else "---"

# Cross-fitting and FSIF rows (7 columns)
cf_row   <- c("No", "No", "$L = 5$", "$L = 5$",
              "No", "$L = 5$", "$L = 5$")
fsif_row <- c("No", "No", "Yes", "Yes",
              "No", "Yes", "Yes")

# Trimming info
trim_row <- c("---", "---",
              "$[0.10, 0.90]$",
              sprintf("%.2f", est$alpha_optimal),
              "---",
              "$[0.10, 0.90]$",
              sprintf("%.2f", est$alpha_optimal))

# Effective sample size after trimming
pG_raw <- est$pG_raw_cf
alpha_opt <- est$alpha_optimal
n_eff_symm <- format(sum(pG_raw >= 0.10 & pG_raw <= 0.90), big.mark = ",")
n_eff_opt  <- format(sum(pG_raw <= (1 - alpha_opt)), big.mark = ",")
n_eff_row <- c(N_str, N_str,
               n_eff_symm, n_eff_opt,
               N_str,
               n_eff_symm, n_eff_opt)

# Covariates row
p_row <- c("---", p_str, p_str, p_str, p_str, p_str, p_str)

# LaTeX table (7 columns)
ncols <- 7

tex <- c(
  sprintf("\\begin{tabular}{l%s}", paste(rep("c", ncols), collapse = "")),
  "\\toprule",
  paste0(" & \\multicolumn{4}{c}{Wald Estimand}",
         " & \\multicolumn{3}{c}{TC Estimand} \\\\"),
  "\\cmidrule(lr){2-5} \\cmidrule(lr){6-8}",
  paste0(" & Wald & Wald$_X$ & DML-Wald & DML-Wald",
         " & TC$_X$ & DML-TC & DML-TC \\\\"),
  paste0(" & & & (symmetric) & (optimal)",
         " & & (symmetric) & (optimal) \\\\"),
  "\\midrule",
  sprintf("Point estimate & %s \\\\",
          paste(est_row, collapse = " & ")),
  sprintf("SE (clustered) & %s \\\\",
          paste(se_clust_row, collapse = " & ")),
  sprintf("95\\%% CI & %s \\\\",
          paste(ci_row, collapse = " & ")),
  "\\midrule",
  sprintf("$DID_D$ (first stage) & \\multicolumn{%d}{l}{%s (primary school completion rate)} \\\\",
          ncols, did_D_str),
  sprintf("$DID_Y$ (reduced form) & \\multicolumn{%d}{l}{%s} \\\\",
          ncols, did_Y_str),
  sprintf("Observations & \\multicolumn{%d}{c}{%s} \\\\",
          ncols, N_str),
  sprintf("Effective $N$ & %s \\\\",
          paste(n_eff_row, collapse = " & ")),
  sprintf("Covariates ($p$) & %s \\\\",
          paste(p_row, collapse = " & ")),
  sprintf("Cross-fitting & %s \\\\",
          paste(cf_row, collapse = " & ")),
  sprintf("FSIF correction & %s \\\\",
          paste(fsif_row, collapse = " & ")),
  sprintf("Trimming ($\\hat{\\alpha}$) & %s \\\\",
          paste(trim_row, collapse = " & ")),
  "\\bottomrule",
  "\\end{tabular}"
)

outfile1 <- file.path(DIR_TABLES, "empirical_duflo.tex")
writeLines(tex, outfile1)
message(sprintf("  Saved: %s", outfile1))

# =============================================================================
# 4. TABLE 2: ESTIMATOR COMPARISON / DECOMPOSITION
# =============================================================================

message("Generating empirical_comparison.tex...")

# Compute differences to highlight what drives variation
wald_vs_waldx <- waldx$estimate - wald$estimate
waldx_vs_drwald <- drw_o$estimate - waldx$estimate
tcx_vs_drtc <- drt_o$estimate - tcx$estimate
wald_vs_tc <- tcx$estimate - wald$estimate
drwald_vs_drtc <- drt_o$estimate - drw_o$estimate

# Confidence interval widths
ci_widths <- sapply(results, function(r) r$ci_upper - r$ci_lower)

tex2 <- c(
  "\\small",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Comparison & $\\hat{\\Delta}_A$ & $\\hat{\\Delta}_B$ & Difference & Interpretation \\\\",
  "\\midrule",
  sprintf("Wald vs.\\ Wald$_X$ & %s & %s & %s & Covariate adjustment \\\\",
          fmt_est(wald$estimate), fmt_est(waldx$estimate), fmt_est(wald_vs_waldx)),
  sprintf("Wald$_X$ vs.\\ DR-Wald & %s & %s & %s & Debiasing (FSIF + CF) \\\\",
          fmt_est(waldx$estimate), fmt_est(drw_o$estimate), fmt_est(waldx_vs_drwald)),
  sprintf("TC$_X$ vs.\\ DR-TC & %s & %s & %s & Debiasing (FSIF + CF) \\\\",
          fmt_est(tcx$estimate), fmt_est(drt_o$estimate), fmt_est(tcx_vs_drtc)),
  sprintf("Wald vs.\\ TC$_X$ & %s & %s & %s & Estimand (Wald vs.\\ TC) \\\\",
          fmt_est(wald$estimate), fmt_est(tcx$estimate), fmt_est(wald_vs_tc)),
  sprintf("DR-Wald vs.\\ DR-TC & %s & %s & %s & Estimand (Wald vs.\\ TC) \\\\",
          fmt_est(drw_o$estimate), fmt_est(drt_o$estimate), fmt_est(drwald_vs_drtc)),
  "\\midrule",
  "\\multicolumn{5}{l}{\\textit{95\\% CI widths}} \\\\",
  "\\midrule",
  sprintf("Wald & \\multicolumn{4}{c}{%.3f} \\\\", ci_widths[1]),
  sprintf("Wald$_X$ & \\multicolumn{4}{c}{%.3f} \\\\", ci_widths[2]),
  sprintf("TC$_X$ & \\multicolumn{4}{c}{%.3f} \\\\", ci_widths[3]),
  sprintf("DR-Wald (optimal) & \\multicolumn{4}{c}{%.3f} \\\\", ci_widths[4]),
  sprintf("DR-TC (optimal) & \\multicolumn{4}{c}{%.3f} \\\\", ci_widths[7]),
  "\\bottomrule",
  "\\end{tabular}"
)

outfile2 <- file.path(DIR_TABLES, "empirical_comparison.tex")
writeLines(tex2, outfile2)
message(sprintf("  Saved: %s", outfile2))

# =============================================================================
# 5. TABLE 3: PRE-TREND / PLACEBO TESTS
# =============================================================================

message("Generating empirical_placebo.tex...")

placebo_emp <- est$placebo_employment
placebo_ps  <- est$placebo_propensity

if (!is.null(placebo_emp) && !is.null(placebo_ps)) {
  tex3 <- c(
    "\\small",
    "\\begin{tabular}{lccc}",
    "\\toprule",
    "Test & Estimate & SE & $p$-value \\\\",
    "\\midrule",
    sprintf("Placebo outcome (employment) & %s & %s & %.3f \\\\",
            add_stars(placebo_emp$estimate, placebo_emp$se),
            fmt_est(placebo_emp$se),
            placebo_emp$pvalue),
    sprintf("\\quad $DID_Y^{\\text{emp}}$ & \\multicolumn{3}{c}{%s} \\\\",
            fmt_est(placebo_emp$did_empY)),
    "\\midrule",
    sprintf("Propensity-split falsification & %s & %s & %.3f \\\\",
            add_stars(placebo_ps$estimate, placebo_ps$se),
            fmt_est(placebo_ps$se),
            placebo_ps$pvalue),
    sprintf("\\quad Median $\\hat{e}(X)$ & \\multicolumn{3}{c}{%.3f} \\\\",
            placebo_ps$pG_median),
    "\\bottomrule",
    "\\end{tabular}"
  )

  outfile3 <- file.path(DIR_TABLES, "empirical_placebo.tex")
  writeLines(tex3, outfile3)
  message(sprintf("  Saved: %s", outfile3))
} else {
  message("  Placebo results not found in estimates. Skipping table.")
}

# =============================================================================
# 6. CONSOLE SUMMARY
# =============================================================================

message("\n========== TABLE SUMMARY ==========")
nms <- c("Wald", "Wald_X", "TC_X",
         "DR-Wald(opt)", "DR-Wald(no-trim)", "DR-Wald(fixed)",
         "DR-TC(opt)", "DR-TC(no-trim)", "DR-TC(fixed)")
for (i in seq_along(results)) {
  r <- results[[i]]
  clust_str <- if (!is.null(r$se_clustered) && !is.na(r$se_clustered)) {
    sprintf("  ClustSE=%s", fmt_est(r$se_clustered))
  } else ""
  message(sprintf("  %18s: %s  SE=%s%s  CI=%s",
                  nms[i], add_stars(r$estimate, r$se),
                  fmt_est(r$se), clust_str, fmt_ci(r$ci_lower, r$ci_upper)))
}

if (!is.null(placebo_emp)) {
  message(sprintf("\n  Placebo (employment): %.4f (p=%.3f)",
                  placebo_emp$estimate, placebo_emp$pvalue))
}
if (!is.null(placebo_ps)) {
  message(sprintf("  Placebo (propensity split): %.4f (p=%.3f)",
                  placebo_ps$estimate, placebo_ps$pvalue))
}

message("\n=== Stage 4 complete ===")
