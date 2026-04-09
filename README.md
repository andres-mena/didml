# didml

**Difference-in-Differences with Machine Learning**

Unified R package for DID estimation with machine learning nuisance estimation. Supports sharp DID, fuzzy DID, and extensible multi-period designs.

## Available Estimators

| Design | IV (Fuzzy) | DML | Estimator | Reference |
|--------|-----------|-----|-----------|-----------|
| 2x2 | No | No | DR-DID | Sant'Anna & Zhao (2020) |
| 2x2 | No | Yes | DML-DID | Chang (2020) |
| 2x2 | Yes | Yes | DML-Wald, DML-TC | Mena (2026) |
| Multi-period | * | * | *Coming soon* | |

## Installation

```r
# install.packages("remotes")
remotes::install_github("andres-mena/didml")

# For sharp DID without DML (optional):
install.packages("DRDID")
```

## Quick Start

```r
library(didml)

# Sharp DID with ML (Chang 2020)
fit_sharp <- didml(Y, D, G, Ti, X, iv = FALSE, dml = TRUE)

# Sharp DID without ML (Sant'Anna & Zhao 2020)
fit_drdid <- didml(Y, D, G, Ti, X, iv = FALSE, dml = FALSE)

# Fuzzy DID with ML (Mena 2026)
fit_fuzzy <- didml(Y, D, G, Ti, X, iv = TRUE, dml = TRUE)

# Print, summarize, plot
print(fit_fuzzy)
summary(fit_fuzzy)
plot(fit_fuzzy)
```

## Bundled Application

```r
# INPRES school construction (Duflo 2001, CD'H 2018)
data(duflo)
fit <- didml(duflo$Y, duflo$D, duflo$G, duflo$Ti, duflo$X,
             iv = TRUE, dml = TRUE, method = "lasso", K = 5,
             cluster = duflo$cluster)
print(fit)
```

## Features

- **Three estimation paths**: Sharp (DRDID), Sharp DML (Chang), Fuzzy DML (Wald/TC)
- **Single entry point**: `didml()` routes to the right estimator via `iv` and `dml` flags
- **Flexible ML**: Lasso, Random Forest, neural networks, OLS, or custom functions
- **Cross-fitting**: K-fold sample splitting for valid post-selection inference
- **Optimal trimming**: Data-driven one-sided propensity score trimming
- **Inference**: Analytical, cluster-robust, or multiplier bootstrap SEs
- **Bundled data**: INPRES school construction application

## Modular API

```r
# Step-by-step for advanced users (fuzzy DID)
nuis <- didml_nuisance(Y, D, G, Ti, X, method = "rf", K = 5, iv = TRUE)
trim <- didml_trim(nuis$pG_raw, method = "auto")
wald <- didml_wald(W, nuis)
inf  <- didml_inference(wald, cluster = district_id, se_type = "cluster")

# Step-by-step (sharp DID with DML)
nuis <- didml_nuisance(Y, D, G, Ti, X, method = "lasso", K = 5, iv = FALSE)
chang <- didml_chang(W, nuis)
inf   <- didml_inference(chang)
```

## Related Software

| Package | Language | Scope |
|---------|----------|-------|
| **didml** (this) | R | DID-specific: sharp, fuzzy, DML, trimming |
| [DRDID](https://github.com/pedrohcgs/DRDID) | R | Sharp DR-DID (Sant'Anna & Zhao 2020) |
| [ddml](https://github.com/thomaswiemann/ddml) | R | General DDML (partial linear, interactive, IV) |
| [ddml](https://github.com/aahrens1/ddml) | Stata | General DDML + stacking (Ahrens et al. 2024) |
| [csdid](https://github.com/d2cml-ai/csdid) | R/Stata | Multi-period DID (Callaway & Sant'Anna 2021) |

## References

- Ahrens, A., Hansen, C.B., Schaffer, M.E. and Wiemann, T. (2024). ddml: Double/debiased machine learning in Stata. *The Stata Journal*, 24(1), 3-45.
- Chang, N.-C. (2020). Double/debiased machine learning for difference-in-differences models. *The Econometrics Journal*, 23(2), 177-191.
- Chernozhukov, V. et al. (2018). Double/Debiased Machine Learning. *Econometrics Journal*, 21(1).
- De Chaisemartin, C. and D'Haultfoeuille, X. (2018). Fuzzy Differences-in-Differences. *Review of Economic Studies*, 85(2).
- Mena, A. (2026). Double Debiased Machine Learning for DID under Imperfect Compliance. Brown University.
- Sant'Anna, P.H.C. and Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics*, 219(1), 101-122.

## License

MIT
