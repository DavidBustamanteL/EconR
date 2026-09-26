# econr 0.1.2

- Added `hatvalues.econ_lm()` so sandwich HC2/HC3 can obtain leverage values directly from `econ_lm` objects.
- Weighted leverage uses the actual WLS design `sqrt(w) * X`.
- Added HC3 compatibility regression tests.

# econr 0.1.1

- `sandwich` is now an imported dependency.
- HC0-HC3 covariance matrices are delegated to `sandwich::vcovHC()`.
- One-way clustered covariance is delegated to `sandwich::vcovCL()`.
- Added `bread.econ_lm()` and `estfun.econ_lm()` so `econ_lm` objects participate directly in the sandwich ecosystem.
- Added equivalence tests against `sandwich` applied to base `lm` objects.

# econr 0.1.0

- Initial Phase-1 implementation of OLS/WLS and model-specific NA handling.
