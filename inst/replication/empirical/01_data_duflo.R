# =============================================================================
# 01_data_duflo.R — Load SUPAS, build CD'H groups, construct quartile dictionary
#
# Replicates CD'H (2018) Section 5.2 group classification on the Duflo (2001)
# INPRES data. Builds a quartile-bin covariate dictionary (p ~ 197) from
# district-level variables. Saves in the same format as 01_data.R so that
# 02_nuisance.R and 03_estimators.R can run unchanged.
#
# Input:  data/raw/inpresdata.dta (Roodman repo)
# Output: data/clean/duflo_analysis.rds
#
# Author: Andres Mena (Brown University)
# =============================================================================

source("code/R/sim/00_config.R")

library(haven)

set.seed(EMPIRICAL_SEED)

# =============================================================================
# 1. LOAD SUPAS DATA
# =============================================================================

message("=== 01_data_duflo.R: Load SUPAS + CD'H group classification ===\n")

raw <- as.data.frame(read_dta("data/raw/inpresdata.dta"))
message(sprintf("Raw: %s obs, %d vars", format(nrow(raw), big.mark = ","), ncol(raw)))

# =============================================================================
# 2. SAMPLE RESTRICTIONS
# =============================================================================

# Both sexes, positive wages, cohort 0 (born 1957-62) or cohort 1 (born 1968-72)
d <- raw[!is.na(raw$lwage) & (raw$treat1b == 1 | raw$treat2b == 1), ]
message(sprintf("Both sexes, wages, cohorts 0+1: %s", format(nrow(d), big.mark = ",")))

d$Ti <- d$treat2b                        # young cohort (born >= 1968)
d$D  <- as.numeric(d$yeduc >= 6)         # completed primary school (binary)
d$Y  <- d$lwage                          # log monthly wages

# =============================================================================
# 3. CD'H SECTION 5.2 GROUP CLASSIFICATION
# =============================================================================

# Bin education into 5 categories (CD'H footnote 12)
d$educ_cat <- cut(d$yeduc, breaks = c(-Inf, 5, 8, 11, 14, Inf),
                  labels = 1:5, right = TRUE)

districts <- unique(d$birthpl)
district_info <- data.frame(birthpl = districts, G_star = NA_integer_,
                            pvalue = NA_real_, delta_D = NA_real_)

for (i in seq_along(districts)) {
  dd <- d[d$birthpl == districts[i], ]
  dd_old <- dd[dd$Ti == 0, ]; dd_young <- dd[dd$Ti == 1, ]

  if (nrow(dd_old) < 5 || nrow(dd_young) < 5) {
    district_info$pvalue[i] <- 1; district_info$delta_D[i] <- 0; next
  }

  tab <- table(factor(dd$educ_cat, levels = 1:5), dd$Ti)
  test <- tryCatch(chisq.test(tab), error = function(e) list(p.value = 1))
  district_info$pvalue[i] <- test$p.value
  district_info$delta_D[i] <- mean(dd_young$D) - mean(dd_old$D)
}

district_info$G_star <- ifelse(
  is.na(district_info$pvalue) | district_info$pvalue > 0.5, 0L,
  ifelse(district_info$delta_D > 0, 1L, -1L)
)

message(sprintf("\nDistrict classification:"))
message(sprintf("  G*=1  (increased): %d", sum(district_info$G_star == 1)))
message(sprintf("  G*=0  (control):   %d", sum(district_info$G_star == 0)))
message(sprintf("  G*=-1 (decreased): %d", sum(district_info$G_star == -1)))

# Merge and keep G*=1 vs G*=0
d <- merge(d, district_info[, c("birthpl", "G_star")], by = "birthpl", all.x = TRUE)
d <- d[d$G_star %in% c(0, 1), ]
d$G <- as.numeric(d$G_star == 1)
rownames(d) <- NULL

N <- nrow(d)
message(sprintf("\nAnalysis sample (G*=1 vs G*=0): N = %s", format(N, big.mark = ",")))
print(table(G = d$G, Ti = d$Ti))

# First stage check
did_D <- mean(d$D[d$G==1 & d$Ti==1]) - mean(d$D[d$G==1 & d$Ti==0]) -
        (mean(d$D[d$G==0 & d$Ti==1]) - mean(d$D[d$G==0 & d$Ti==0]))
did_Y <- mean(d$Y[d$G==1 & d$Ti==1]) - mean(d$Y[d$G==1 & d$Ti==0]) -
        (mean(d$Y[d$G==0 & d$Ti==1]) - mean(d$Y[d$G==0 & d$Ti==0]))
message(sprintf("\nDID_D: %.4f", did_D))
message(sprintf("DID_Y: %.4f", did_Y))
message(sprintf("Wald:  %.3f (CD'H Table 3: 0.140)", did_Y / did_D))

# =============================================================================
# 4. BUILD QUARTILE DICTIONARY
# =============================================================================

cont_vars <- c("nin", "ch71", "en71", "wsppc", "dens71", "moldyed")
n_bins <- 4

bin_list <- list()
bin_var_labels <- list()

for (v in cont_vars) {
  x <- d[[v]]
  breaks <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(breaks) < 3) breaks <- sort(unique(c(min(x, na.rm=T)-1, unique(x[!is.na(x)]), max(x, na.rm=T)+1)))
  bins <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bins[is.na(bins)] <- 1L
  n_actual <- length(unique(bins))
  bm <- matrix(0, nrow = nrow(d), ncol = n_actual)
  for (b in 1:n_actual) bm[, b] <- as.numeric(bins == b)
  colnames(bm) <- paste0(v, "_Q", 1:n_actual)
  bm <- bm[, -1, drop = FALSE]
  bin_list[[v]] <- bm
  bin_var_labels[[v]] <- colnames(bm)
  message(sprintf("  %s: %d bins → %d indicators", v, n_actual, ncol(bm)))
}

X_bins <- do.call(cbind, bin_list)
X_java  <- matrix(d$java, ncol = 1, dimnames = list(NULL, "java"))
X_urban <- matrix(d$urban, ncol = 1, dimnames = list(NULL, "urban"))
X_atomic <- cbind(X_java, X_urban, X_bins)
atom_names <- colnames(X_atomic)
n_atoms <- ncol(X_atomic)

# Pairwise interactions (exclude within-variable bin products)
int_cols <- list()
for (i in 1:(n_atoms - 1)) {
  for (j in (i + 1):n_atoms) {
    same_var <- FALSE
    for (v in cont_vars) {
      if (atom_names[i] %in% bin_var_labels[[v]] &&
          atom_names[j] %in% bin_var_labels[[v]]) { same_var <- TRUE; break }
    }
    if (same_var) next
    ix <- X_atomic[, i] * X_atomic[, j]
    if (var(ix, na.rm = TRUE) > 1e-10)
      int_cols[[paste0(atom_names[i], ":", atom_names[j])]] <- ix
  }
}
X_int <- do.call(cbind, int_cols)

# Full dictionary
X_matrix <- cbind(d[, cont_vars, drop = FALSE], X_atomic, X_int)
for (j in 1:ncol(X_matrix)) {
  nas <- is.na(X_matrix[, j])
  if (any(nas)) X_matrix[nas, j] <- median(X_matrix[!nas, j])
}
col_vars <- apply(X_matrix, 2, var, na.rm = TRUE)
X_matrix <- as.matrix(X_matrix[, col_vars > 1e-10, drop = FALSE])

variable_names <- colnames(X_matrix)
message(sprintf("\nDictionary: p = %d", ncol(X_matrix)))

# =============================================================================
# 5. SAVE (same format as 01_data.R for pipeline compatibility)
# =============================================================================

dir.create(DIR_EMPIRICAL, recursive = TRUE, showWarnings = FALSE)

out <- list(
  Y              = d$Y,
  D              = d$D,
  G              = d$G,
  Ti             = d$Ti,
  STATEFIP       = d$birthpl,  # district ID for clustering
  X_matrix       = X_matrix,
  variable_names = variable_names,
  df             = d,
  metadata       = list(
    source         = "SUPAS 1995 (Duflo 2001), CD'H groups (Section 5.2)",
    N              = N,
    p              = ncol(X_matrix),
    districts      = length(unique(d$birthpl)),
    group_counts   = table(district_info$G_star),
    did_D          = did_D,
    did_Y          = did_Y,
    wald_naive     = did_Y / did_D,
    cohort_old     = "born 1957-1962 (treat1b)",
    cohort_young   = "born 1968-1972 (treat2b)",
    D_definition   = "completed primary school (yeduc >= 6)",
    Y_definition   = "log monthly wages (lwage)"
  )
)

outfile <- file.path(DIR_EMPIRICAL, "duflo_analysis.rds")
saveRDS(out, outfile)
message(sprintf("\nSaved: %s (%.1f MB)", outfile,
                file.size(outfile) / 1024^2))

message("\n=== Stage 1 (Duflo) complete ===")
