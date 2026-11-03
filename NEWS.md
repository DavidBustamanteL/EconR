
# EconR 0.1.5

- Added `econ_multinomial()` for Multinomial Logit and Multinomial Probit regression.
- Implemented Multinomial Logit estimation using `nnet::multinom()` as the estimation backend.
- Implemented Multinomial Probit estimation with fixed or estimated error covariance matrices, using multivariate normal probabilities from `mvtnorm`.
- Added classical, HC0, HC1, CR0, CR1, CR1p, CR1S, and Stata-style covariance estimators.
- Added analytical observation-level score contributions through `estfun.econ_multinomial()`.
- Added support for observation weights, subsetting, missing-value handling, missing-weight warnings, and cluster alignment.
- Added parameter inference, confidence intervals, joint Wald tests, and model summary statistics.
- Implemented in-sample and out-of-sample probability predictions, predicted classes, and prediction standard errors.
- Added average marginal effects for continuous and categorical predictors, including weighted estimation and uncertainty quantification.
- Integrated multinomial models with standard R methods, including `coef()`, `vcov()`, `summary()`, `predict()`, `fitted()`, `logLik()`, `nobs()`, and `confint()`.
- Added compatibility with `broom::tidy()`, `broom::glance()`, and `broom::augment()`.
- Validated Multinomial Logit estimates against `nnet` and Multinomial Probit estimates against `mlogit`.
- Added independent numerical checks for covariance estimation, score contributions, prediction probabilities, and marginal effects.
- Completed 485 passing assertions in the multinomial test suite, with zero failures, warnings, or skipped expectations.

# EconR 0.1.4

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
