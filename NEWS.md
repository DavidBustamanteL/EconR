# EconR 0.1.3

- Added `econ_ordered()` for Ordered Logit and Ordered Probit regression, using `MASS::polr()` as the estimation backend.
- Implemented classical, HC0, HC1, CR0, and CR1 covariance estimators for ordered models.
- Added analytical observation-level score contributions through `estfun.econ_ordered()`.
- Added support for observation weights, subsetting, missing-value handling, and cluster alignment.
- Integrated ordered models with standard R methods, including `coef()`, `vcov()`, `summary()`, `predict()`, and `logLik()`.
- Added compatibility with `broom::tidy()`, `broom::glance()`, `broom::augment()`, and `modelsummary`.
- Added 155 automated assertions covering estimation, robust and clustered covariance, analytical scores, predictions, reporting, and input validation.
- Validated ordered-model estimates and covariance matrices against `MASS` and `sandwich`.

# EconR 0.1.2

- Added `econ_binary()` for binary Logit and Probit regression.
- Added `hatvalues.econ_lm()` so sandwich HC2/HC3 can obtain leverage values directly from `econ_lm` objects.
- Weighted leverage uses the actual WLS design `sqrt(w) * X`.
- Added HC3 compatibility regression tests.

# EconR 0.1.1

- `sandwich` is now an imported dependency.
- HC0–HC3 covariance matrices are delegated to `sandwich::vcovHC()`.
- One-way clustered covariance is delegated to `sandwich::vcovCL()`.
- Added `bread.econ_lm()` and `estfun.econ_lm()` so `econ_lm` objects participate directly in the sandwich ecosystem.
- Added equivalence tests against `sandwich` applied to base `lm` objects.

# EconR 0.1.0

- Initial Phase-1 implementation of OLS/WLS and model-specific NA handling.
