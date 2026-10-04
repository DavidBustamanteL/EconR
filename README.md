
# EconR

`EconR` is a developing R package for econometric estimation with a consistent interface and model-specific missing-data handling.

The package aims to provide a unified framework for econometric models while preserving familiar R formula syntax and compatibility with established R infrastructure.

## Current Version

**EconR 0.1.2**

The current release implements linear and binary-response models:

- Ordinary Least Squares (OLS)
- Weighted Least Squares (WLS)
- Binary Logit
- Binary Probit

Both model families support robust and clustered covariance estimation, prediction, inference, and standard R model methods.

## Features

### Linear Models

The `econ_lm()` function supports:

- Ordinary Least Squares (OLS)
- Weighted Least Squares (WLS)
- R formula syntax
- Factors and custom contrasts
- Transformations and interactions
- Offsets
- Subsetting
- Rank-deficient models
- Zero-weight observations

### Binary Response Models

The `econ_binary()` function supports:

- Binary Logit (`link = "logit"`)
- Binary Probit (`link = "probit"`)
- Weighted maximum-likelihood estimation
- R formula syntax
- Factors, transformations, and interactions
- Subsetting and missing-data handling
- Classical, heteroskedasticity-robust, and clustered covariance estimation
- In-sample and out-of-sample prediction
- Predicted probabilities and linear predictors
- Log-likelihood, AIC, and BIC
- Coefficient inference using z-statistics and normal confidence intervals

The dependent variable must be numeric and binary, containing only `0` and `1`.

### Missing-Data Handling

`EconR` constructs the estimation sample from the variables actually required by the specified model.

Missing values in unrelated columns therefore do not remove observations from the estimation sample.

Supported `NA` behavior includes:

- `na.omit`
- `na.exclude`
- `na.fail`

Missing values in model inputs such as weights and cluster identifiers are incorporated into the same estimation-sample logic.

This approach is shared by `econ_lm()` and `econ_binary()`.

### Covariance Estimation

Both `econ_lm()` and `econ_binary()` support:

- `"classical"`
- `"HC0"`
- `"HC1"`
- `"HC2"`
- `"HC3"`
- `"CR0"`
- `"CR2"`
- `"stata"`

Robust covariance estimation builds on established infrastructure from `sandwich` and `clubSandwich`.

`econ_lm` and `econ_binary` objects implement the estimating-function and bread-matrix interfaces required by `sandwich`.

For example:

```r
fit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "HC3"
)

vcov(fit)
```

### Standard R Methods

Both linear and binary-response models support:

- `coef()`
- `vcov()`
- `predict()`
- `summary()`
- `model.matrix()`
- `nobs()`
- `hatvalues()`

Additional linear-model methods include:

- `confint()`
- `anova()`
- `fitted()`
- `residuals()`
- `weights()`
- `formula()`

Additional binary-response functionality includes:

- `logLik()`
- `AIC()`
- `BIC()`

Compatibility methods are also provided for:

- `sandwich`
- `broom`
- `modelsummary`

Both model families support `broom::tidy()` and `broom::glance()`, allowing fitted models to be incorporated into standard reporting workflows.

## Installation

The development version can be installed directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("DavidBustamanteL/EconR")
```

To install the tagged `v0.1.2` release specifically, once published:

```r
remotes::install_github(
    "DavidBustamanteL/EconR@v0.1.2"
)
```

## Examples

### 1. Ordinary Least Squares (OLS)

```r
library(EconR)

fit = econ_lm(
    mpg ~ wt * factor(am) + hp,
    data = mtcars,
    vcov = "HC3"
)

summary(fit)
confint(fit)
anova(fit)

predict(
    fit,
    newdata = mtcars[1:3, ]
)
```

### 2. Clustered Linear Model

Clustered covariance estimation can be requested directly:

```r
fit_cluster = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    vcov = "CR2",
    cluster = cyl
)

summary(fit_cluster)
```

### 3. Binary Logit

```r
fit_logit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "HC3"
)

summary(fit_logit)

# Predicted probabilities
predict(
    fit_logit,
    type = "response"
)

# Predictions for new observations
predict(
    fit_logit,
    newdata = mtcars[1:3, ],
    type = "response"
)
```

### 4. Binary Probit

The same interface can be used for Probit estimation:

```r
fit_probit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "probit",
    vcov = "HC1"
)

summary(fit_probit)

predict(
    fit_probit,
    type = "response"
)
```

### 5. Clustered Binary Logit

```r
fit_logit_cluster = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "CR2",
    cluster = cyl
)

summary(fit_logit_cluster)
```

### 6. Weighted Binary Logit

Observation weights can be supplied directly:

```r
mtcars$w = mtcars$disp / mean(mtcars$disp)

fit_weighted = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    weights = w,
    link = "logit",
    vcov = "HC1"
)

summary(fit_weighted)
```

### 7. Compatibility with broom and modelsummary

```r
library(broom)
library(modelsummary)

# Coefficient table
tidy(fit_logit)

# Model-level statistics
glance(fit_logit)

# Regression table
modelsummary(
    list(
        "Logit" = fit_logit,
        "Probit" = fit_probit
    ),
    stars = TRUE,
    output = "default"
)
```

### 8. Direct Sandwich Compatibility

The fitted model can also be used directly with `sandwich`:

```r
sandwich::vcovHC(
    fit_logit,
    type = "HC3"
)
```

## Validation

Version `0.1.2` has been tested against established R implementations, including `lm()`, `glm()`, `estimatr`, `sandwich`, and `clubSandwich`.

Validation covers coefficient estimation, covariance matrices, predictions, missing-data handling, inference, and compatibility with standard R model-reporting tools.

The automated test suite currently contains:

| Model family | Passing tests |
|---|---:|
| Linear models (OLS/WLS) | 403 |
| Binary models (Logit/Probit) | 139 |
| **Total** | **542** |

The complete package successfully passes `R CMD check`:

```text
0 errors | 0 warnings | 0 notes
```

## Design Principle

The central design principle of `EconR` is that the estimation sample should be determined by the model being estimated.

An observation is excluded because information required by that model is missing—not because an unrelated variable elsewhere in the data contains a missing value.

`EconR` therefore owns model construction, estimation-sample handling, and estimation semantics, while established packages such as `sandwich` and `clubSandwich` provide specialized covariance-estimation infrastructure.

The same design philosophy is applied consistently across linear and binary-response models.

## Development Status

`EconR` is under active development.

Version `0.1.2` establishes the validated OLS/WLS and Logit/Probit frameworks, including robust and clustered inference, model-specific missing-data handling, and compatibility with standard R reporting tools.

Future development will extend this interface to additional econometric model classes.

## Author

David Bustamante Lazo

## Repository

https://github.com/DavidBustamanteL/EconR
