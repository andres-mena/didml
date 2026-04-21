# didml 0.4.0

- Add full paper replication package under `inst/replication/`.
  - `sim/` — Monte Carlo pipeline (5 stages) producing every simulation
    table and figure in the paper (main text and appendix).
  - `empirical/` — INPRES / Duflo (2001) application producing every
    empirical table and figure, including the faithful port of de
    Chaisemartin & D'Haultfoeuille (2018) Table 3 via Roodman's `late()`
    Mata function.
  - `shared/` — canonical estimator implementations used by both pipelines.
  - `paper_tables/`, `paper_figures/` — frozen copies of every exhibit as
    it appears in the paper PDF, for byte-level diff against
    user-regenerated output.
  - `README.md` — maps every paper exhibit to its producer script.

# didml 0.3.1

- Fix two critical bugs flagged by CRAN prechecks.
- Add validation guards on user-supplied learners.

# didml 0.3.0

- Stacking of multiple ML learners.
- Repeated cross-fitting.
- Separate learner selection for each nuisance component.
- Comprehensive test suite (143 tests).
- Simulation DGP helper `dgp_didml()` matching paper Section 5.

# didml 0.2.x

- DML-Wald and DML-TC estimators for fuzzy DID.
- Optimal propensity trimming rule.
- Clustered standard errors.
- Bundled INPRES data.
