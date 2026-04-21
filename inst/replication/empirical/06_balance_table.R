# =============================================================================
# 06_balance_table.R -- Descriptive statistics by (G,T) cell
#
# Generates a publication-quality LaTeX table with means and SDs
# of key variables for each (G,T) cell.
#
# Input:  data/clean/duflo_analysis.rds
# Output: output/tables/empirical_balance.tex
#
# Author:  Andres Mena (Brown University)
# Project: fuzzyDDML
# =============================================================================

source("code/R/sim/00_config.R")

message("=== 06_balance_table.R: Descriptive statistics ===\n")

dat <- readRDS("data/clean/duflo_analysis.rds")
df  <- dat$df

dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)

# Variables and labels
vars   <- c("lwage", "D", "yeduc", "nin", "ch71", "en71",
            "wsppc", "dens71", "moldyed", "java", "urban")
labels <- c("Log monthly wages ($Y$)", "Primary completion ($D$)",
            "Years of education", "INPRES intensity",
            "Children pop.\\ 1971", "Enrollment rate 1971",
            "Water/sanitation spending", "Pop.\\ density 1971",
            "Avg.\\ education (old cohort)", "Java", "Urban")

# Cell subsets
cells <- list(
  c00 = df[df$G == 0 & df$Ti == 0, ],
  c01 = df[df$G == 0 & df$Ti == 1, ],
  c10 = df[df$G == 1 & df$Ti == 0, ],
  c11 = df[df$G == 1 & df$Ti == 1, ]
)
cell_n <- sapply(cells, nrow)

# Compute means and SDs
make_row <- function(v, lab) {
  means <- sapply(cells, function(d) mean(d[[v]], na.rm = TRUE))
  sds   <- sapply(cells, function(d) sd(d[[v]], na.rm = TRUE))
  # Format: 2 decimal places for most, 3 for small values
  fmt <- ifelse(abs(means) < 1, "%.3f", "%.2f")
  mean_str <- mapply(sprintf, fmt, means)
  sd_str   <- mapply(sprintf, fmt, sds)
  row_mean <- paste0(lab, " & ", paste(mean_str, collapse = " & "), " \\\\")
  row_sd   <- paste0(" & ", paste(paste0("(", sd_str, ")"), collapse = " & "), " \\\\")
  c(row_mean, row_sd)
}

# Build table
lines <- character()
lines <- c(lines, "\\begin{tabular}{lcccc}")
lines <- c(lines, "\\toprule")
lines <- c(lines, " & \\multicolumn{2}{c}{Control ($G=0$)} & \\multicolumn{2}{c}{Treatment ($G=1$)} \\\\")
lines <- c(lines, "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}")
lines <- c(lines, sprintf(" & $T=0$ & $T=1$ & $T=0$ & $T=1$ \\\\"))
lines <- c(lines, sprintf(" & Old cohort & Young cohort & Old cohort & Young cohort \\\\"))
lines <- c(lines, "\\midrule")

# Outcome and treatment
lines <- c(lines, "\\multicolumn{5}{l}{\\textit{Panel A: Outcome and treatment}} \\\\")
lines <- c(lines, "\\midrule")
for (i in 1:3) {
  lines <- c(lines, make_row(vars[i], labels[i]))
}

lines <- c(lines, "\\\\[0.3em]")
lines <- c(lines, "\\multicolumn{5}{l}{\\textit{Panel B: District characteristics}} \\\\")
lines <- c(lines, "\\midrule")
for (i in 4:length(vars)) {
  lines <- c(lines, make_row(vars[i], labels[i]))
}

lines <- c(lines, "\\midrule")
n_str <- formatC(cell_n, format = "d", big.mark = "{,}")
lines <- c(lines, paste0("$N$ & ", paste(n_str, collapse = " & "), " \\\\"))
lines <- c(lines, "\\bottomrule")
lines <- c(lines, "\\end{tabular}")

tex <- paste(lines, collapse = "\n")
writeLines(tex, file.path(DIR_TABLES, "empirical_balance.tex"))
message(sprintf("  Saved: %s", file.path(DIR_TABLES, "empirical_balance.tex")))

message("\n=== Stage 6 complete ===")
