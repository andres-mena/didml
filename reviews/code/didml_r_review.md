# didml R Package -- Code Quality Review

**Reviewer:** Senior R Developer (automated review)
**Date:** 2026-04-09
**Package version:** 0.2.0
**Files reviewed:** 11 files in `R/`

---

## Executive Summary

The `didml` package is well-structured and demonstrates strong software engineering practices for an econometrics research package. The code is cleanly organized, thoroughly documented with roxygen2, and handles edge cases carefully. Key strengths include disciplined RNG state preservation, informative error messages, and a clean separation of concerns across files. The main issues are: one potentially serious domain bug (variable naming vs. semantics in `didml_chang()`), a `message()` call in library code outside verbose mode, and several medium-priority items around robustness and consistency.

**Overall package score: 82/100**

---

## File-by-File Review

---

### 1. `didml.R` -- Main Entry Point

**Score: 88/100**

This is the user-facing orchestrator. Excellent roxygen2 documentation, clear three-path branching (DRDID / Chang / Wald-TC), and good input validation delegation.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 116 | Low | Function Design | `match.arg(design, ...)` is called but `design` default is `"2x2"`, which is fine. However, `estimator` on line 250 only validates for the fuzzy path. If a user passes `estimator = "wald"` with `iv = FALSE, dml = TRUE`, they reach the Chang path silently. No error or warning. | Add explicit validation: if `!iv && dml && estimator != "auto" && estimator != "chang"`) warn or error that the estimator was overridden. |
| 173 | Low | Script Structure | Class vector `c("didml", "didml_sharp")` puts the generic class first and the specific subclass second. R dispatches on the *first* class in the vector. | Swap to `c("didml_sharp", "didml")` so that any future `print.didml_sharp` method gets dispatched first. Same issue on lines 245 and 333. |
| 262 | Low | Domain Correctness | `DID_D_hat` is computed for trimming but `pT` on line 268 is recomputed as `mean(Ti)` while `nuis` already contains `pT` from the nuisance step. These should agree but it is redundant. | Use `nuis$pT` or remove the redundant assignment. |
| 141-145 | Low | Function Design | The DRDID path passes `D` to `.didml_drdid()`, but DRDID ignores `D` (it uses `G` as its treatment). The `D` argument is accepted but unused in the wrapper. Not a bug per se but could confuse readers. | Add a brief comment explaining that `D` is passed for the compliance check only. |

**Strengths:**
- Clear step-by-step verbose messages (gated behind `verbose`).
- Consistent return structure across all three paths.
- `stringsAsFactors = FALSE` used consistently in `data.frame()`.
- Clean `match.call()` capture.

---

### 2. `nuisance.R` -- Cross-Fitted Nuisance Estimation

**Score: 83/100**

The workhorse of the DML pipeline. Handles both sharp and fuzzy paths with clean fold-based cross-fitting. The NA-not-impute discipline is excellent.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 102 | Medium | Domain Correctness | Hard clipping of propensity to `[0.001, 0.999]` happens *before* the trimming step in `didml.R`. This silently modifies extreme propensities before the user's chosen trimming rule takes effect. The trimming function never sees the true extreme values. | Document this pre-clipping clearly in the roxygen, or defer all clipping to `didml_trim()`. |
| 186 | Medium | Domain Correctness | Same hard-clipping issue as line 102, in the fuzzy path. | Same fix. |
| 149-153 | Low | Domain Correctness | The cell `"11"` (G=1, Ti=1) is intentionally excluded from conditional expectation estimation (nuisance is not needed for the treated-in-post cell). This is correct per the theory but nowhere documented in comments. | Add a brief comment: `# Cell (1,1) excluded: nuisance not needed for treated-post obs`. |
| 240 | Medium | Reproducibility | Inside the TC-specific loop, `L_ind` is constructed using `seq_len(N)` on every fold iteration. This is correct but wasteful: `L_ind` is the same for all folds of a given cell. | Move `L_ind` computation outside the fold loop. |
| 240 | Medium | Domain Correctness | The propensity `pi_{dgt}` is trained on the full sample (`folds != k`), not within the control group. `L_ind` is `as.numeric(seq_len(N) %in% idx)` where `idx` filters on D, G, Ti. This means the propensity model is `P(D=d, G=g, T=t | X)` estimated from all observations. This is correct for a joint propensity, but the comment on line 239 says "Propensity for this (d,g,t) cell" which could be clearer. | Clarify comment to "Joint cell-membership propensity P(D=d, G=g, T=t | X) trained on full sample." |
| 87 | Low | Comment Quality | `pT <- mean(Ti)` -- no comment explaining that this relies on the assumption T independent of (G, X). The @details roxygen mentions it but an inline note would help. | Add: `# pT = marginal P(T=1), valid under T indep. of (G, X)`. |

**Strengths:**
- Explicit NA initialization vectors (never rely on partial assignment).
- Clean separation of sharp vs. fuzzy nuisance requirements.
- Warning when NAs are detected in core predictions.
- `pG_raw` naming convention clearly distinguishes from trimmed values.

---

### 3. `scores.R` -- DML-Wald, DML-TC, DML-Chang Score Functions

**Score: 79/100**

Contains the three core estimator functions. The Wald and TC implementations are careful and well-documented. However, `didml_chang()` has a naming issue that could cause confusion, and there is one potential domain-correctness concern.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 232 | **Critical** | Domain Correctness | `p_0 <- mean(G)` is labeled `# Pr(G = 1) = P(D=1) in Chang's notation`. But `p_0` is used as the *denominator* divisor. In Chang (2020), the weight denominator includes `Pr(D=0)` (probability of being untreated in Chang's notation), which equals `1 - mean(G)` in our notation. The variable is named `p_0` suggesting "probability of control group" but it is assigned `mean(G)` which is the probability of the *treated* group. This is internally consistent (the formula on line 241 uses `p_0` as `Pr(G=1)` in the denominator), but the variable name `p_0` is misleading -- it reads as "prob of G=0." | Rename `p_0` to `p_G` or `pG_bar` to avoid ambiguity. The formula itself appears correct: `denom_w = pT * (1 - pT) * mean(G) * (1 - e(X))` matches Chang's Eq 3.2 where the weight normalizes by `P(D=1) * (1 - e(X))` (in Chang's D = our G). |
| 234 | Low | Error Handling | The guard `p_0 < 1e-10 || p_0 > 1 - 1e-10` checks for degenerate group composition. The threshold `1e-10` is extremely tight -- with N=1000, having just 1 treated observation gives `p_0 = 0.001` which passes this check. | Consider a more practical threshold, e.g., `sum(G) < 10 || sum(1-G) < 10`. |
| 260-261 | Medium | Domain Correctness | `ATT <- mean(psi_i[ok])` averages over non-NA scores. If some scores are NA (trimmed or failed nuisance), the effective N changes but the estimate is still an average over surviving observations. This is correct if NAs are MCAR relative to the estimand, but could introduce bias if NAs are systematically from extreme propensities. | Add a warning when the fraction of NA scores exceeds, say, 10%, similar to the nuisance NA check. |
| 67-70 | Low | Comment Quality | The Wald score formula is well-commented at the top but the individual lines of the score computation lack inline annotation mapping to paper equations. | Add equation references inline, e.g., `# Eq. (3.1) in Mena (2026)`. |
| 126-129 | Medium | Domain Correctness | In `didml_tc()`, the propensities `pred$pi_101` etc. are floored with `pmax(..., 0.01)`. This is a different floor than the `0.001` used in nuisance.R line 245. Inconsistent trimming floors across the pipeline. | Unify the minimum propensity floor to a single constant, e.g., `.MIN_PROPENSITY <- 0.001`. |
| 147-151 | Low | Domain Correctness | The TC numerator formula is complex. No reference to the specific equation in the paper is given. | Add a reference comment: `# TC numerator: Theorem 2 of Mena (2026)`. |

**Strengths:**
- Weak first-stage detection with informative warnings.
- Non-finite weight guards with NA (not zero) replacement.
- `wald_naive` diagnostic returned alongside debiased estimate.
- Clean Wald ratio structure: `theta_hat = sum(num_i) / sum(den_i)`.

---

### 4. `inference.R` -- Standard Errors and Confidence Intervals

**Score: 85/100**

Clean implementation of three SE methods. The bootstrap uses score-based multiplier (Chernozhukov et al. 2018), which is correct and efficient.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 44-47 | **High** | Console Output Hygiene | `message("Note: \`cluster\` provided; switching to cluster-robust SEs.")` is called unconditionally when `cluster != NULL && se_type != "cluster"`. This prints to the console in library code outside any `verbose` guard. | Wrap in a `verbose` parameter (requires adding `verbose` to function signature) or use `warning()` instead of `message()`. |
| 58 | Medium | Domain Correctness | Cluster-robust variance: `V_hat <- (S / (S - 1)) * sum(cluster_sums^2) / (N^2 * G_hat^2)`. The finite-sample correction `S/(S-1)` is standard. However, for ratio estimands (Wald, TC), `G_hat` should be the Jacobian, and the formula assumes the scores `psi_i` are already centered at the estimated theta. This is correct for this package's score definitions. | No fix needed, but add a comment: `# Scores already centered at theta_hat; G_hat = Jacobian`. |
| 82-83 | Low | Reproducibility | The bootstrap loop `for (b in seq_len(B))` is sequential. For large B, this is slow. | Consider documenting that `B` should be moderate (e.g., 999 or 1999) or offering parallelism in a future version. |
| 89-90 | Medium | Domain Correctness | Analytical variance: `se_hat <- sqrt(V_hat / N)`. For the Chang estimator, `V_hat` from `didml_chang()` is `mean((psi_i - ATT)^2)`, which is already divided by N in the mean. So `se_hat = sqrt(V_hat / N)` gives the correct `sqrt(E[psi^2] / N)`. For Wald/TC, `V_hat = mean(psi^2) / G_hat^2`, same logic. Consistent and correct. | No fix needed. |
| 66-74 | Low | Reproducibility | RNG state preservation is correctly implemented with `on.exit()`. Good practice. | None -- this is well done. |

**Strengths:**
- Three SE methods cleanly separated.
- RNG state preserved and restored via `on.exit()`.
- Rademacher multiplier bootstrap is the correct modern approach.
- Clean conflict resolution (cluster provided but se_type mismatch).

---

### 5. `trim.R` -- Propensity Score Trimming

**Score: 84/100**

Implements the data-driven trimming rule from the paper. Clean grid-search structure.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 35-36 | Low | Domain Correctness | For `method = "fixed"`, trimming is symmetric: `pmax(alpha, pmin(1-alpha, pG))`. But the data-driven method (line 62-67) is one-sided (only caps at `1-alpha`). The inconsistency between fixed and auto trimming behavior is not documented. | Document in roxygen that `fixed` applies symmetric trimming while `auto` applies one-sided upper trimming per Corollary 2.1. |
| 44-45 | Low | Comment Quality | `# Without first-stage predictions, use Corollary 2.2 (constant g)` -- good reference but which paper? | Add: `# Corollary 2.2 of Mena (2026)`. |
| 62-63 | Low | Domain Correctness | The one-sided trimming selects only on `pG_raw <= (1 - a)`. It does not also trim observations with very low `pG_raw` (near 0). In theory, low propensity (overlap violations from the other direction) could also inflate variance. | Consider whether lower-bound trimming should also be searched in `auto` mode, or document the one-sided choice. |
| 64 | Low | Error Handling | `if (sum(sel) < 10)` uses the same magic number as `.MIN_CELL_SIZE`. | Replace with `.MIN_CELL_SIZE` for consistency. |
| 47-49 | Low | Error Handling | `tryCatch` around loess is good, but the fallback (`g_hat = rep(mean(DID_D), N)`) does not warn the user that loess failed. | Add a message or warning when falling back to constant `g_hat`. |

**Strengths:**
- Clean three-mode interface (none / fixed / auto).
- Grid search is simple and transparent.
- Returns `keep` indicator for downstream diagnostics.
- Loess failure is caught gracefully.

---

### 6. `methods.R` -- S3 Methods (print, summary, plot, confint, tidy)

**Score: 87/100**

Comprehensive set of S3 methods. Clean `ggplot2` usage with Brown University color scheme (#012169).

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 104 | Low | Console Output Hygiene | `print(x$estimates, row.names = FALSE, digits = 4)` uses base `print.data.frame()` which outputs to the console directly. This is fine in a `print.summary` method, but the `digits` argument to `print.data.frame` controls significant digits, not decimal places. | Consider using `format()` or `sprintf()` for more control, matching the style in `print.didml()`. |
| 142-154 | Medium | Function Design | `tidy.didml()` is exported but is not registered as an S3 method for `generics::tidy()` or `broom::tidy()`. The NAMESPACE shows `export(tidy.didml)` but no `S3method(tidy, didml)`. If `broom` or `generics` is loaded, `tidy(fit)` will not dispatch to this method. | Add `@method tidy didml` to the roxygen block and ensure `S3method(tidy, didml)` appears in NAMESPACE. Or import `generics::tidy` and register properly. |
| 114-133 | Low | Function Design | `confint.didml()` recomputes CIs using `z * se` which may differ from stored `ci_lower`/`ci_upper` if a different `level` is requested. This is correct behavior but could surprise users if they expect `confint(fit)` to match `fit$estimates$ci_lower`. | Add a note in documentation that `confint()` recomputes intervals at the requested level. |
| 181-201 | Low | Figure Quality | The propensity plot uses `bins = 40` which may be too many for small samples. | Consider using `ggplot2::geom_density()` as an alternative or adaptive bin count. |
| 168 | Low | Function Design | `plot.didml()` only accepts `type = "propensity"` and `type = "scores"`. No `type = "coefficient"` or similar for DRDID path. The DRDID path has no nuisance and will error on the propensity plot. | Either add a guard for the DRDID case (already done on line 175-177, good), or add a `type = "estimates"` option for all paths. |

**Strengths:**
- Clean `ggplot2` usage with `theme_minimal()`.
- Expression-based axis labels (proper math rendering).
- Trimming threshold shown on propensity plot with annotation.
- `invisible(x)` returned from print methods.
- `confint()` properly handles `parm` subsetting.

---

### 7. `drdid_wrap.R` -- DRDID Wrapper

**Score: 82/100**

Clean wrapper around `DRDID::drdid_rc()`. The D-G mapping is well documented.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 41 | Medium | Domain Correctness | `X_mat <- cbind(1, as.matrix(X))` adds an intercept column. However, `DRDID::drdid_rc()` with its default settings may already add an intercept internally, leading to collinearity. | Check DRDID documentation -- if `drdid_rc()` adds its own intercept, remove the `cbind(1, ...)`. If not, this is correct. |
| 33-39 | Low | Error Handling | Compliance rate warning threshold is 0.95 (5% non-compliance). This is somewhat arbitrary. | Consider making this threshold configurable or documenting the choice. |
| 48-49 | Low | Function Design | `nboot = if (use_boot) B else NULL` -- when `use_boot` is FALSE, passing `nboot = NULL` to DRDID may trigger a warning depending on the DRDID version. | Pass `nboot = 0L` or omit entirely using `do.call()`. |
| 25 | Low | Function Design | The `cluster` parameter is accepted but never used. The comment says "reserved for future use." | Either implement cluster support or remove the parameter to avoid confusion. |
| 60-61 | Low | Polish | `if (!is.null(drdid_out$att.inf.func)) drdid_out$att.inf.func else NULL` is equivalent to `drdid_out$att.inf.func` since accessing a NULL list element returns NULL. | Simplify to `scores = drdid_out$att.inf.func`. |

**Strengths:**
- Clear mapping documentation: DRDID's `D` = our `G`, DRDID's `post` = our `Ti`.
- Compliance rate diagnostic is a thoughtful addition.
- Dependency check with `requireNamespace()`.

---

### 8. `compat.R` -- Deprecated Function Wrappers

**Score: 90/100**

Simple, clean deprecation layer. Uses `.Deprecated()` correctly.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 16-22 | Low | Function Design | `fddml()` maps `estimand` to `estimator` but doesn't handle `estimand = "chang"` or `estimand = "drdid"` which are valid for `didml()` but not the original `fddml()`. No error, just potentially confusing. | Add validation: `estimand <- match.arg(estimand, c("wald", "tc", "both"))` before dispatching. |
| 26-29 | Low | Function Design | `fddml_nuisance()` passes all args via `...` but doesn't explicitly set `iv = TRUE` (the original default). If the new `didml_nuisance()` changes defaults, behavior could silently change. | Explicitly set `iv = TRUE` in the call. |

**Strengths:**
- Clean one-file deprecation strategy.
- Uses `.Deprecated()` with proper package attribution.
- All deprecated functions delegate to new equivalents.

---

### 9. `utils.R` -- Internal Helpers

**Score: 86/100**

Contains validation, fold creation, and the ML fitting dispatch. This is the most complex internal file.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 58-59 | Low | Console Output Hygiene | `warning("Cell has < ", .MIN_CELL_SIZE, ...)` fires inside a loop that runs K * num_cells times. For K=5 and 3 cells, this could produce 15 warnings if many cells are small. | Consider collecting warnings and emitting a single summary warning after all folds are processed (requires refactoring the caller in `nuisance.R`). |
| 81-82 | Low | Reproducibility | `glmnet::cv.glmnet(..., nfolds = 5)` uses internal cross-validation with no seed control. The glmnet fold assignment is random and not governed by the outer `seed` argument. | Set a local seed before `cv.glmnet()` calls or use the `foldid` argument. |
| 94-96 | Low | Function Design | Random Forest uses hardcoded `ntree = 500L`, `mtry = max(1, floor(p/3))`, `nodesize = 5L`. These are reasonable defaults but not tunable. | Consider accepting `...` for method-specific arguments, or document these defaults in the roxygen for `didml_nuisance()`. |
| 110-111 | Low | Function Design | Neural network uses hardcoded `size = 5L`, `decay = 0.01`, `maxit = 500L`. Same tunability concern as RF. | Same suggestion as above. |
| 4 | Low | Polish | `.MIN_CELL_SIZE <- 10L` is defined at file scope. Good practice, but it's only referenced by name inside `utils.R` itself. The magic number `10` also appears in `trim.R` line 64 and `nuisance.R` (implicitly via the function call). | Use `.MIN_CELL_SIZE` in `trim.R` as well for consistency. |
| 16 | Low | Error Handling | `if (!all(G %in% 0:1))` checks values but not type. If G is character "0"/"1", it would pass (since `"0" %in% 0:1` is FALSE due to type coercion). Actually this would correctly error. No issue. | None needed. |
| 24-37 | Low | Reproducibility | `.make_folds()` preserves and restores RNG state. Excellent. | None -- well implemented. |

**Strengths:**
- `.MIN_CELL_SIZE` is a named constant, not a magic number.
- `tryCatch` wraps all ML estimation with informative error propagation.
- `.validate_inputs()` checks lengths, binary constraints, and NAs comprehensively.
- RNG state preservation in `.make_folds()` is correct.
- Family-aware estimation (binomial vs. gaussian) across all methods.

---

### 10. `data.R` -- Dataset Documentation

**Score: 92/100**

Thorough roxygen2 documentation for the `duflo` dataset.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 9 | Low | Comment Quality | `X` is described as "Numeric matrix (N x 197)" but the exact construction (quartile-binned district characteristics plus pairwise interactions) could benefit from a reference to the specific replication code or paper section. | Add: "See Appendix B of Mena (2026) for construction details." |
| 41-44 | Low | Polish | The example uses `method = "lasso"` and `K = 5` which may be slow for a quick demo. | Consider adding `\dontrun{}` or reducing K for the example. |

**Strengths:**
- Clear description of sample construction and cohort definitions.
- Source and references properly cited.
- `@format` describes all list components.
- Group classification methodology documented.

---

### 11. `globals.R` -- Global Variable Declarations

**Score: 95/100**

Standard `utils::globalVariables()` to suppress R CMD check NOTEs for NSE.

| Line | Severity | Category | Issue | Suggested Fix |
|------|----------|----------|-------|---------------|
| 1-9 | Low | Polish | Some variables listed may not actually trigger NOTEs (e.g., `fold`, `component`, `value`, `n_label`). These might be from earlier iterations of the code. | Run `R CMD check` and remove any variables that don't produce NOTEs. |

**Strengths:**
- Correct use of `globalVariables()` for ggplot2 NSE variables.
- All relevant variables appear to be covered.

---

## Cross-File Issues

| Severity | Category | Issue | Files Affected | Suggested Fix |
|----------|----------|-------|----------------|---------------|
| **Critical** | Domain Correctness | Variable `p_0` in `didml_chang()` (scores.R:232) is assigned `mean(G)` but named to suggest `Pr(G=0)`. While the formula is correct, the misleading name risks bugs in future modifications. | scores.R | Rename to `pG_bar` or `p_1`. |
| **High** | Console Output Hygiene | `message()` on inference.R:46 fires unconditionally in library code, violating the principle that packages should be silent unless verbose mode is on. | inference.R | Gate behind `verbose` or convert to `warning()`. |
| Medium | Domain Correctness | Propensity floor is `0.001` in nuisance.R but `0.01` in scores.R (TC weights). Inconsistent. | nuisance.R, scores.R | Define a single constant `.MIN_PROPENSITY` in utils.R and use throughout. |
| Medium | Domain Correctness | Pre-clipping of propensity to `[0.001, 0.999]` in nuisance.R happens before `didml_trim()`. The trimming function never sees the true extremes. | nuisance.R, trim.R, didml.R | Either defer all clipping to `didml_trim()` or document the two-stage process. |
| Medium | Function Design | `tidy.didml()` is exported but not registered as S3 method. Will not dispatch via `tidy()` generic. | methods.R, NAMESPACE | Add `@method tidy didml` and ensure S3 registration. |
| Low | Polish | Class vectors are ordered `c("didml", "didml_*")` but R dispatches on the first element. | didml.R | Swap to `c("didml_sharp", "didml")` or `c("didml_fuzzy", "didml")`. |
| Low | Reproducibility | `glmnet::cv.glmnet()` internal folds are not seeded, introducing minor non-reproducibility in nuisance estimates even when `seed` is set. | utils.R | Use `foldid` argument or set a local seed. |

---

## Score Summary

| File | Score | Key Issue |
|------|-------|-----------|
| `didml.R` | 88 | Minor: estimator validation gap, class order |
| `nuisance.R` | 83 | Medium: pre-clipping before trimming, L_ind inefficiency |
| `scores.R` | 79 | Critical: misleading `p_0` name; Medium: inconsistent propensity floors |
| `inference.R` | 85 | High: unconditional `message()` in library code |
| `trim.R` | 84 | Low: asymmetric trim docs, magic number |
| `methods.R` | 87 | Medium: `tidy.didml` S3 registration |
| `drdid_wrap.R` | 82 | Medium: possible double intercept |
| `compat.R` | 90 | Low: missing arg validation in deprecated wrappers |
| `utils.R` | 86 | Low: unseeded glmnet CV, hardcoded ML hyperparameters |
| `data.R` | 92 | Low: example speed |
| `globals.R` | 95 | Low: possibly stale variable names |

**Weighted average: 83/100**

---

## Priority Action Items

### Must Fix (Critical/High)

1. **scores.R:232** -- Rename `p_0` to `pG_bar` or `p_treated` to avoid a variable name that reads as "probability of G=0" but actually stores `Pr(G=1)`. The formula is correct, but the naming is a bug waiting to happen.

2. **inference.R:44-47** -- The `message()` call fires unconditionally when `cluster` is provided but `se_type != "cluster"`. In library code, this should be gated behind verbose or converted to a `warning()`.

### Should Fix (Medium)

3. **nuisance.R + scores.R** -- Unify propensity floor constants (`0.001` vs `0.01`) into a single `.MIN_PROPENSITY` constant in `utils.R`.

4. **nuisance.R:102,186** -- Document or reconsider pre-clipping of propensity scores before the trimming step.

5. **methods.R:142** -- Register `tidy.didml()` as a proper S3 method for the `generics::tidy` generic.

6. **drdid_wrap.R:41** -- Verify whether `DRDID::drdid_rc()` adds its own intercept; if so, remove `cbind(1, ...)`.

### Nice to Have (Low)

7. Swap class vector ordering in `didml.R`.
8. Seed `glmnet::cv.glmnet()` internal folds in `utils.R`.
9. Use `.MIN_CELL_SIZE` consistently in `trim.R`.
10. Add equation references in score computation code.

---

*Report generated 2026-04-09. No source files were modified.*
