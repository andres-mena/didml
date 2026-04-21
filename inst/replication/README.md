# Replication Package — Mena (2026), "Debiased Machine Learning for Fuzzy DID"

This directory contains the scripts that produce every table and figure in the
paper's simulation (Section 5) and empirical (Section 6) content.

## Directory layout

```
inst/replication/
├── README.md             # this file
├── sim/                  # Monte Carlo simulations (Section 5 + Appendix)
├── empirical/            # INPRES application (Section 6 + Appendix)
├── shared/               # shared estimator implementations used by both pipelines
├── paper_tables/         # frozen .tex files as they appear in the paper
└── paper_figures/        # frozen .pdf files as they appear in the paper
```

## Paper output → producer script

| Paper exhibit | File | Producer script |
|---|---|---|
| Table — main simulation results (Table A.2) | `sim_main.tex` | `sim/04_tables.R` |
| Table — ML-learner comparison (Table A.3) | `sim_ml_comparison.tex` | `sim/07_ml_comparison_table.R` |
| Table — trimming-rule comparison (Table A.4) | `sim_trim_rule_summary.tex` | `sim/trim_rule_summary.R` |
| Figure — DML-Wald trimming bias-variance (Fig. 1) | `sim_trimming_wald.pdf` | `sim/06_trimming_analysis.R` |
| Table — INPRES returns to schooling (Table 1) | `empirical_duflo.tex` | `empirical/04_tables.R` |
| Table — descriptive/balance statistics | `empirical_balance.tex` | `empirical/06_balance_table.R` |
| Table — four-learner robustness (Table B.3) | `empirical_4learner.tex` | `empirical/four_learner_table.R` |
| Table — classifier comparison vs CDH, Duflo | `empirical_classifier.tex` | static (hand-written) |
| Figure — propensity under Lasso (Fig. B.1a) | `empirical_propensity_lasso.pdf` | `empirical/07_figures_appendix.R` |
| Figure — propensity under Ridge (Fig. B.1b) | `empirical_propensity_ridge.pdf` | `empirical/07_figures_appendix.R` |
| Figure — sensitivity to $\alpha$, primary (Fig. B.2a) | `empirical_sensitivity_primary.pdf` | `empirical/08_figures_sensitivity.R` |
| Figure — sensitivity to $\alpha$, HS (Fig. B.2b) | `empirical_sensitivity_highschool.pdf` | `empirical/08_figures_sensitivity.R` |

## Pipelines

### Simulations (`sim/`)

Five-stage pipeline. Each stage reads the output of the previous and writes
its own. Stages 01 and 02 are expensive; 03–07 are cheap.

```
00_config.R                  # shared paths, seeds, hyperparameters (single source of truth)
01_dgp.R                     # Stage 1: draw M replications per sample size (EXPENSIVE)
02_nuisance.R                # Stage 2: cross-fitted ML nuisances (EXPENSIVE)
03_estimators.R              # Stage 3: seven estimators per replication
04_tables.R                  # Stage 4: main-text LaTeX tables
05_figures.R                 # Stage 5: main-text figures
06_trimming_analysis.R       # bias-variance curve for the trimming threshold
07_ml_comparison_table.R     # Lasso vs RF vs NN robustness
trim_rule_summary.R          # symmetric vs data-driven trimming-rule comparison
```

Run (from repo root, with a sourced `00_config.R` defining paths):

```bash
Rscript sim/01_dgp.R
Rscript sim/02_nuisance.R --method lasso
Rscript sim/02_nuisance.R --method rf
Rscript sim/02_nuisance.R --method nn
Rscript sim/03_estimators.R
Rscript sim/04_tables.R
Rscript sim/05_figures.R
Rscript sim/06_trimming_analysis.R
Rscript sim/07_ml_comparison_table.R
Rscript sim/trim_rule_summary.R
```

Stages 01 and 02 take multiple hours on a workstation. Stages 03–07 complete
in seconds to minutes.

### Empirical (`empirical/`)

Uses the 1995 Indonesian SUPAS microdata from Duflo (2001), already bundled in
the `didml` package as `data(duflo)`. Replication starts from
`data/raw/inpresdata.dta` (the Roodman / Duflo extract), not from the raw
SUPAS files.

```
01_data_duflo.R                    # build analysis sample, CDH classifier G*
02_nuisance.R                      # cross-fitted ML first step
03_estimators.R                    # DML-Wald, DML-TC, plug-in variants
04_tables.R                        # Table 1 (main empirical table)
06_balance_table.R                 # descriptive balance table
07_figures_appendix.R              # propensity-density figures
08_figures_sensitivity.R           # sensitivity to trimming threshold
four_learner_table.R               # Lasso / Ridge / Elastic / RF comparison
replication_cdh_table3.R           # faithful port of CDH (2018) Table 3
```

`replication_cdh_table3.R` is a direct port of David Roodman's Stata/Mata
`late()` function from `github.com/droodman/Duflo-2001`. It reproduces CDH
(2018) Table 3 to within one bootstrap standard error on the same SUPAS
extract.

### Shared (`shared/`)

```
04_estimators.R   # compute_dd_wald, compute_dd_tc — used by both pipelines
```

The canonical estimator implementations in the `didml` R package itself
(`R/scores.R`, `R/didml.R`) have been validated against this script.

## Frozen outputs (`paper_tables/`, `paper_figures/`)

Exact copies of the files referenced in the paper PDF, as of the package
version date. These are committed so a user can diff their locally-regenerated
output against the paper version without running the full pipeline.

## Data

The simulation pipeline needs no external data. The empirical pipeline needs
`inpresdata.dta` (Duflo/Roodman SUPAS extract, public):

- Source: `github.com/droodman/Duflo-2001/tree/master/data`.
- Variables used: `yeduc, p105, p504thn, lhwage, birthpl, java, urban, weight`.
- Sample after cohort and non-missing-wage filters: $N = 22{,}414$ observations in $201$ districts of birth.

## Citation

If you use this replication package, please cite:

> Mena, A. (2026). "Debiased Machine Learning for Fuzzy Differences-in-Differences." *Working paper, Brown University*.

and the software:

```
citation("didml")
```
