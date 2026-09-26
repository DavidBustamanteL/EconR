# econr

A developing R econometrics package with a consistent interface and model-specific missing-data handling.

## Phase 1

- OLS and WLS
- R formula syntax, factors, transformations and interactions
- model-specific `NA` omission
- classical covariance plus HC0-HC3 and one-way clustered covariance via `sandwich`
- direct compatibility with `sandwich::vcovHC()` and `sandwich::vcovCL()`
- `coef()`, `vcov()`, `fitted()`, `residuals()`, `predict()`, `summary()`

```r
fit <- econ_lm(mpg ~ wt * factor(am) + hp, mtcars, vcov = "HC3")
summary(fit)

# The econr model can also be passed directly to sandwich:
sandwich::vcovHC(fit, type = "HC3")
```

The design principle is that missing values in unrelated columns do not remove an observation from a model. `econr` owns model construction and estimation-sample handling; `sandwich` supplies established robust covariance infrastructure.
