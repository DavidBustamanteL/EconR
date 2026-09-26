# EconR

`EconR` is a developing R package for econometric estimation with a consistent interface and model-specific missing-data handling.

The package aims to provide a unified framework for econometric models while preserving familiar R formula syntax and compatibility with established R infrastructure.

## Current Version

**EconR 0.1.2**

The current release implements the core linear-model framework, including OLS and WLS estimation, robust and clustered covariance estimators, prediction, inference, and standard model methods.

## Features

### Linear Models

- Ordinary Least Squares (OLS)
- Weighted Least Squares (WLS)
- R formula syntax
- Factors and custom contrasts
- Transformations and interactions
- Offsets
- Subsetting
- Rank-deficient models
- Zero-weight observations

### Missing-Data Handling

`EconR` constructs the estimation sample from the variables actually required by the specified model.

Missing values in unrelated columns therefore do not remove observations from the estimation sample.

Supported `NA` behavior includes:

- `na.omit`
- `na.exclude`
- `na.fail`

Missing values in model inputs such as weights and cluster identifiers are incorporated into the same estimation-sample logic.

### Covariance Estimation

`econ_lm()` currently supports:

- `"classical"`
- `"HC0"`
- `"HC1"`
- `"HC2"`
- `"HC3"`
- `"CR0"`
- `"CR2"`
- `"stata"`

Robust covariance estimation builds on established infrastructure from `sandwich` and `clubSandwich`.

`econ_lm` objects can also be passed directly to functions such as:

```r
sandwich::vcovHC(fit, type = "HC3")
```

### Standard R Methods

`econ_lm` objects support:

- `coef()`
- `vcov()`
- `confint()`
- `anova()`
- `fitted()`
- `residuals()`
- `predict()`
- `summary()`
- `model.matrix()`
- `weights()`
- `nobs()`
- `hatvalues()`
- `formula()`

Compatibility methods are also provided for:

- `sandwich`
- `broom`
- `modelsummary`

## Installation

The development version can be installed directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("DavidBustamanteL/EconR")
```

To install the validated `v0.1.2` release specifically:

```r
remotes::install_github(
    "DavidBustamanteL/EconR@v0.1.2"
)
```

## Example

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

The fitted model can also be used directly with `sandwich`:

```r
sandwich::vcovHC(
    fit,
    type = "HC3"
)
```

## Design Principle

The central design principle of `EconR` is that the estimation sample should be determined by the model being estimated.

An observation is excluded because information required by that model is missing—not because an unrelated variable elsewhere in the data contains a missing value.

`EconR` therefore owns model construction, estimation-sample handling, and estimation semantics, while established packages such as `sandwich` and `clubSandwich` provide specialized covariance-estimation infrastructure.

## Development Status

`EconR` is under active development.

Version `0.1.2` establishes the validated OLS/WLS core. Future development will extend the same interface and estimation-sample semantics to additional econometric model classes.

## Author

David Bustamante Lazo

## Repository

https://github.com/DavidBustamanteL/EconR