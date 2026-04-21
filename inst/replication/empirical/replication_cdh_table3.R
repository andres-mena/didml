# =============================================================================
# replication_cdh_table3.R — Replicate de Chaisemartin & D'Haultfoeuille (2018)
# RES Section 5.2, Table 3 on Duflo/INPRES — FAITHFUL port of the Mata late()
# function from Roodman's archive (Code/de Chaisemartin and d'Haultfoeuille
# 2017 simulation.do, lines 28–89).
#
# CRITICAL: CDH's Wald-DID is a POOLED estimator over G ∈ {−1, 0, +1}, not a
# two-group contrast on G ∈ {0, +1}. From late() line 58:
#
#   Wald_x[1] = n_{G=1}·(ΔEY1 − ΔEY0) + n_{G=-1}·(ΔEY0 − ΔEY_neg1)
#
# divided by (ΔΔED + ΔΔED_neg1), where
#   ΔΔED       = n_{G=1}·(ΔED1 − ΔED0)
#   ΔΔED_neg1  = n_{G=-1}·(ΔED0 − ΔED_neg1)
#
# D is ordered (yeduc, 0–18 after recoding 19→18). Y = lhwage.
#
# Target: Wald-DID = 0.140 (SE 0.015).
#
# Input:  data/raw/inpresdata.dta
# Output: reviews/data/2026-04-20_cdh_replication_results.md
# =============================================================================

suppressMessages(library(haven))

cat("=== CDH (2018) Table 3 faithful replication ===\n")
cat("D = yeduc (ordered 0-18), Y = lhwage, Wald-DID pools G ∈ {-1, 0, +1}.\n\n")

d <- as.data.frame(read_dta("data/raw/inpresdata.dta"))
d$yeduc[d$yeduc == 19] <- 18                       # Roodman line 130

d$age74 <- 74 - d$p504thn
d$young <- as.integer(d$age74 <= 6)
d$old   <- as.integer(d$age74 >= 12 & d$age74 <= 17)
s <- d[(d$young == 1 | d$old == 1) & !is.na(d$lhwage), ]
cat(sprintf("Sample after (young|old) & !is.na(lhwage): N = %d\n", nrow(s)))

# -------- CDH classifier ----------
# Per Roodman lines 95-123: chi-sq on full yeduc x young, signed by delta(mean yeduc)
districts <- sort(unique(s$birthpl))
info <- data.frame(birthpl = districts, G_star = 0L, pvalue = NA_real_, delta_yeduc = NA_real_)
for (i in seq_along(districts)) {
  dd <- s[s$birthpl == districts[i], ]
  dd_o <- dd[dd$young == 0, ]; dd_y <- dd[dd$young == 1, ]
  if (nrow(dd_o) == 0 || nrow(dd_y) == 0) { info$G_star[i] <- NA_integer_; next }
  tab <- table(dd$yeduc, dd$young)
  info$pvalue[i]      <- tryCatch(suppressWarnings(chisq.test(tab)$p.value), error = function(e) 1)
  info$delta_yeduc[i] <- mean(dd_y$yeduc) - mean(dd_o$yeduc)
  if (!is.na(info$pvalue[i]) && info$pvalue[i] < 0.5) {
    info$G_star[i] <- ifelse(info$delta_yeduc[i] > 0,  1L,
                       ifelse(info$delta_yeduc[i] < 0, -1L, 0L))
  }
}

cat(sprintf("\nDistrict classification (chi-sq thresh 0.5, sign by delta mean yeduc):\n"))
cat(sprintf("  G* = -1: %d  [CDH: 97]\n",  sum(info$G_star == -1, na.rm = TRUE)))
cat(sprintf("  G* =  0: %d  [CDH: 64]\n",  sum(info$G_star ==  0, na.rm = TRUE)))
cat(sprintf("  G* = +1: %d  [CDH: 123]\n", sum(info$G_star ==  1, na.rm = TRUE)))

s <- merge(s, info[, c("birthpl", "G_star")], by = "birthpl")
# Keep all G* in {-1, 0, +1} (CDH never drop G=-1)
s <- s[!is.na(s$G_star), ]
cat(sprintf("\nFull analysis sample (G* in {-1,0,+1}): N = %d  [CDH: 30,828]\n", nrow(s)))

# -------- late() Mata port ----------
# D = yeduc, Y = lhwage, T = young, G = G_star
D <- s$yeduc; Y <- s$lhwage; T <- s$young; G <- s$G_star

mask <- function(g, t) G == g & T == t

# Line 46-47: ΔED by group
DED0     <- mean(D[mask( 0, 1)]) - mean(D[mask( 0, 0)])
DED1     <- mean(D[mask( 1, 1)]) - mean(D[mask( 1, 0)])
DED_neg1 <- mean(D[mask(-1, 1)]) - mean(D[mask(-1, 0)])

# Line 50-54: ΔEY by group (and EY11, EY_neg11)
DEY0     <-  mean(Y[mask( 0, 1)]) - mean(Y[mask( 0, 0)])
EY11     <-  mean(Y[mask( 1, 1)])
DEY1     <-  EY11 - mean(Y[mask( 1, 0)])
EY_neg11 <-  mean(Y[mask(-1, 1)])
DEY_neg1 <-  EY_neg11 - mean(Y[mask(-1, 0)])

# Line 56-57: weighted ΔΔED
n1   <- sum(G ==  1)
n_n1 <- sum(G == -1)
DDED      <- n1   * (DED1 - DED0)
DDED_neg1 <- n_n1 * (DED0 - DED_neg1)

# Line 58: Wald-DID numerator (pooled)
numW <- n1 * (DEY1 - DEY0) + n_n1 * (DEY0 - DEY_neg1)
denW <- DDED + DDED_neg1

wald_did <- numW / denW

cat(sprintf("\n--- Wald-DID (late() Mata line 58, pooled over G = +1 and G = -1) ---\n"))
cat(sprintf("  n_{G=+1} = %d\n", n1))
cat(sprintf("  n_{G=-1} = %d\n", n_n1))
cat(sprintf("  ΔED1 = %.4f, ΔED0 = %.4f, ΔED_neg1 = %.4f\n", DED1, DED0, DED_neg1))
cat(sprintf("  ΔEY1 = %.4f, ΔEY0 = %.4f, ΔEY_neg1 = %.4f\n", DEY1, DEY0, DEY_neg1))
cat(sprintf("  ΔΔED = %.4f, ΔΔED_neg1 = %.4f\n", DDED, DDED_neg1))
cat(sprintf("  numerator (Wald_x[1]) = %.4f\n", numW))
cat(sprintf("  denominator           = %.4f\n", denW))
cat(sprintf("\nWald-DID = %.4f   [CDH Table 3: 0.140, SE 0.015]\n", wald_did))

# -------- Wald-TC (late() Mata lines 60-86) ----------
# dY0[i] = E[Y | D=d_i, G=0, T=1] - E[Y | D=d_i, G=0, T=0] for each D level
# sum_Q[i] = mean over units in G=1, T=0 at D=d_i of Q(Y00_d, Y01_d, x)
# Wald-TC pools G=+1 and G=-1 similarly via line 86.

minY <- min(Y)

# F2^-1 (F1(x)) — the Q function from Mata lines 14-26
Q_fn <- function(y1, y2, x, val_inf) {
  if (length(y1) == 0 || length(x) == 0) return(rep(NA_real_, length(x)))
  fdr <- sapply(x, function(xv) mean(y1 <= xv))
  n2  <- ceiling(length(y2) * fdr)
  indic_nul <- as.integer(fdr == 0)
  y2_sorted <- sort(y2)
  # Mata: val_inf * indic_nul + (n2 > 0) * sort(y2)[n2 + indic_nul]
  idx <- n2 + indic_nul
  idx <- pmax(pmin(idx, length(y2_sorted)), 1)
  val_inf * indic_nul + (n2 > 0) * y2_sorted[idx]
}

D_levels <- sort(unique(D))
dY0         <- rep(0, length(D_levels))
sum_Q       <- rep(0, length(D_levels))
sum_Q_minus <- rep(0, length(D_levels))

for (k in seq_along(D_levels)) {
  d_k <- D_levels[k]
  Yd  <- Y[D == d_k]; Gd <- G[D == d_k]; Td <- T[D == d_k]
  Y00_d <- Yd[Gd == 0 & Td == 0]
  Y01_d <- Yd[Gd == 0 & Td == 1]
  dY0[k] <- if (length(Y00_d) > 0 && length(Y01_d) > 0) mean(Y01_d) - mean(Y00_d) else 0
  x_plus  <- Yd[Gd ==  1 & Td == 0]
  x_minus <- Yd[Gd == -1 & Td == 0]
  sum_Q[k]       <- if (length(x_plus)  > 0) mean(Q_fn(Y00_d, Y01_d, x_plus,  minY)) else 0
  sum_Q_minus[k] <- if (length(x_minus) > 0) mean(Q_fn(Y00_d, Y01_d, x_minus, minY)) else 0
}

# Weights by D-level for G=+1, T=0 and G=-1, T=0
wt_pos <- sapply(D_levels, function(d_k) sum(G ==  1 & T == 0 & D == d_k))
wt_neg <- sapply(D_levels, function(d_k) sum(G == -1 & T == 0 & D == d_k))

wmean <- function(x, w) sum(x * w) / sum(w)

# Line 79-80: G=+1 Wald-TC / Wald-CIC (before pooling)
wtc_plus  <- (DEY1 - wmean(dY0, wt_pos)) / DED1
wcic_plus <- (EY11 - wmean(sum_Q, wt_pos)) / DED1

# Line 83-84: G=-1 raw numerators
wtc_neg_num  <- DEY_neg1 - wmean(dY0,         wt_neg)
wcic_neg_num <- EY_neg11 - wmean(sum_Q_minus, wt_neg)

# Line 86: pool — (ΔΔED * Wald_plus + ΔΔED_neg1 * W_neg_num/ΔED_neg1) / (ΔΔED + ΔΔED_neg1)
wald_tc  <- (DDED * wtc_plus  + DDED_neg1 * wtc_neg_num  / DED_neg1) / denW
wald_cic <- (DDED * wcic_plus + DDED_neg1 * wcic_neg_num / DED_neg1) / denW

cat(sprintf("\nWald-TC  = %.4f   [CDH Table 3: 0.101]\n", wald_tc))
cat(sprintf("Wald-CIC = %.4f   [CDH Table 3: 0.099]\n", wald_cic))

# -------- Save report ----------
out_path <- "reviews/data/2026-04-20_cdh_replication_results.md"
lines <- c(
  "# CDH (2018) Table 3 Faithful Replication",
  "",
  sprintf("**Run date:** %s", Sys.Date()),
  sprintf("**Data:** `data/raw/inpresdata.dta`"),
  "**Port of:** `late()` Mata function from Roodman's `Code/de Chaisemartin and d'Haultfoeuille 2017 simulation.do` (lines 28–89).",
  "",
  "## Classifier",
  "",
  sprintf("| G* | ours | CDH |"),
  sprintf("|---|---|---|"),
  sprintf("| -1 | %d | 97 |",  sum(info$G_star == -1, na.rm = TRUE)),
  sprintf("| 0  | %d | 64 |",  sum(info$G_star ==  0, na.rm = TRUE)),
  sprintf("| +1 | %d | 123 |", sum(info$G_star ==  1, na.rm = TRUE)),
  "",
  sprintf("Full analysis sample N = %d (target 30,828).", nrow(s)),
  "",
  "## Results — Table 3",
  "",
  "| Estimator | ours | CDH Table 3 (SE) | gap |",
  "|---|---|---|---|",
  sprintf("| Wald-DID | %.4f | 0.140 (0.015) | %+.4f |", wald_did, wald_did - 0.140),
  sprintf("| Wald-TC  | %.4f | 0.101 (0.017) | %+.4f |", wald_tc,  wald_tc  - 0.101),
  sprintf("| Wald-CIC | %.4f | 0.099 (0.017) | %+.4f |", wald_cic, wald_cic - 0.099),
  "",
  "## Mata code faithfully ported (key lines)",
  "",
  "- Line 31: `// D can be finitely supported and ordered, not only binary.`",
  "- Lines 46-47: `ΔED0 = mean(D, G==0 & T) - mean(D, G==0 & !T)` (and G==1).",
  "- Lines 56-57: `ΔΔED = sum(G==1) * (ΔED1 - ΔED0)`; `ΔΔED_neg1 = sum(G==-1) * (ΔED0 - ΔED_neg1)`.",
  "- Line 58: `Wald_x[1] = sum(G==1) * (ΔEY1 - ΔEY0) + sum(G==-1) * (ΔEY0 - ΔEY_neg1)`.",
  "- Line 88: `return Wald_x / (ΔΔED + ΔΔED_neg1)`.",
  "",
  "## Notes",
  "",
  "- D = `yeduc` (ordered, 19 levels 0–18 after recoding 19→18 per Roodman line 130).",
  "- Y = `lhwage` (log hourly wage).",
  "- Both sexes; no sex filter.",
  "- SE requires cluster bootstrap at `birthpl`, 100 reps, seed 0124810 (not run here).",
  "- CDH Wald-DID is return per additional year of schooling (ordered D)."
)
writeLines(lines, out_path)
cat(sprintf("\nSaved: %s\n", out_path))
