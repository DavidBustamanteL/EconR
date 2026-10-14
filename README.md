
# EconR

`EconR` is a developing R package for econometric estimation with a consistent interface and model-specific missing-data handling.

The package aims to provide a unified framework for econometric models while preserving familiar R formula syntax and compatibility with established R infrastructure.

## Current Version

**EconR 0.1.3**

The current version implements three model families:

- Ordinary Least Squares (OLS)
- Weighted Least Squares (WLS)
- Binary Logit and Probit
- Ordered Logit and Probit

The package supports classical, robust, and clustered covariance estimation, prediction, inference, and integration with standard R model-reporting tools.

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
- Classical, heteroskedasticity-robust, and clustered covariance estimation

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

### Ordered Response Models

The `econ_ordered()` function supports:

- Ordered Logit (`link = "logit"`)
- Ordered Probit (`link = "probit"`)
- Weighted maximum-likelihood estimation
- R formula syntax
- Factors, transformations, and interactions
- Subsetting and missing-data handling
- Classical covariance estimation
- Heteroskedasticity-robust covariance estimation (`HC0`, `HC1`)
- Cluster-robust covariance estimation (`CR0`, `CR1`)
- Coefficient and threshold estimation
- Coefficient inference using z-statistics and normal p-values
- Model deviance and AIC
- Integration with `broom` and standard R prediction methods

The dependent variable must be an ordered factor with at least three levels.

Ordered models are estimated using `MASS::polr()`, while `EconR` provides its own covariance estimation and reporting methods.

### Missing-Data Handling

`EconR` constructs the estimation sample from the variables actually required by the specified model.

Missing values in unrelated columns therefore do not remove observations from the estimation sample.

Supported `NA` behavior includes:

- `na.omit`
- `na.exclude`
- `na.fail`

Missing values in model inputs such as weights and cluster identifiers are handled within the model-specific estimation workflow.

This design principle applies across linear, binary-response, and ordered-response models.

### Covariance Estimation

The available covariance estimators depend on the model family.

| Covariance estimator | `econ_lm()` | `econ_binary()` | `econ_ordered()` |
|---|:---:|:---:|:---:|
| `"classical"` | Yes | Yes | Yes |
| `"HC0"` | Yes | Yes | Yes |
| `"HC1"` | Yes | Yes | Yes |
| `"HC2"` | Yes | Yes | — |
| `"HC3"` | Yes | Yes | — |
| `"HC4"` | Yes | Yes | — |
| `"HC4m"` | Yes | Yes | — |
| `"HC5"` | Yes | Yes | — |
| `"CR0"` | Yes | Yes | Yes |
| `"CR1"` | Yes | Yes | Yes |
| `"CR1p"` | Yes | Yes | — |
| `"CR1S"` | Yes | Yes | — |
| `"CR2"` | Yes | Yes | — |
| `"CR3"` | Yes | Yes | — |
| `"stata"` | Yes | Yes | — |

For linear and binary-response models, EconR supports all seven heteroskedasticity-consistent (HC) covariance estimators available through `sandwich::vcovHC()` and all six cluster-robust (CR) covariance estimators available through `clubSandwich::vcovCR()`.

HC estimators account for heteroskedasticity, while CR estimators additionally allow within-cluster error dependence. The `"stata"` option provides Stata-style cluster-robust covariance estimation.

#### Inference Conventions

For binary Logit and Probit models, the reference distribution depends on the covariance estimator:

- **Classical and HC0–HC5:** Asymptotic normal inference (z-tests).
- **CR0, CR1, CR1p, CR1S, CR3, and stata:** Student's t inference with `G - 1` degrees of freedom, where `G` is the number of clusters.
- **CR2:** Student's t inference with coefficient-specific Satterthwaite degrees of freedom, calculated using `clubSandwich`.

These conventions are applied consistently to coefficient p-values and confidence intervals in `summary()`, `tidy()`, and `confint()`.

The availability of covariance estimators for ordered-response models differs from that for linear and binary-response models.

Robust covariance estimation builds on established infrastructure from `sandwich` and `clubSandwich`.

Linear and binary-response models implement estimating-function and bread-matrix interfaces for compatibility with `sandwich`.

Ordered-response models provide observation-level score contributions and covariance estimation for the supported robust and clustered estimators.

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

### Significance Symbols

Model summaries provide significance symbols based on coefficient p-values.

By default:

| Symbol | Significance level |
|---|---|
| `*` | p < 0.10 |
| `**` | p < 0.05 |
| `***` | p < 0.01 |

An optional additional significance threshold is available through:

```r
summary(fit, extra.significance = TRUE)
```

When enabled, coefficients with p-values below `0.001` receive the `+` symbol.

| Symbol | Significance level |
|---|---|
| `*` | p < 0.10 |
| `**` | p < 0.05 |
| `***` | p < 0.01 |
| `+` | p < 0.001 |

Each coefficient receives the symbol corresponding to its most stringent satisfied threshold.

These symbols affect the printed model summary. Numerical p-values remain available for reporting through `broom` and other tools.

### Standard R Methods

`EconR` provides standard R model interfaces, including:

- `coef()`
- `vcov()`
- `predict()`
- `summary()`
- `model.matrix()`
- `nobs()`

Additional linear-model methods include:

- `confint()`
- `anova()`
- `fitted()`
- `residuals()`
- `weights()`
- `formula()`
- `hatvalues()`

Additional binary-response functionality includes:

- `logLik()`
- `AIC()`
- `BIC()`
- `hatvalues()`

Ordered-response models retain compatibility with relevant methods inherited from `MASS::polr`.

Compatibility is provided for established reporting and covariance-estimation infrastructure, including:

- `sandwich`
- `broom`
- `modelsummary`

The package supports integration with `broom::tidy()` for coefficient reporting, including ordered-response coefficients and thresholds.

## Installation

The development version can be installed directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("DavidBustamanteL/EconR")
```

To install the tagged `v0.1.3` release specifically, once published:

```r
remotes::install_github(
    "DavidBustamanteL/EconR@v0.1.3"
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

### 7. Ordered Logit

The `MASS::housing` dataset provides an example of an ordered dependent variable.

```r
library(MASS)

fit_ordered_logit = econ_ordered(
    Sat ~ Infl + Type + Cont,
    data = housing,
    weights = Freq,
    link = "logit",
    vcov = "HC1"
)

summary(fit_ordered_logit)
```

The model estimates regression coefficients and thresholds separating the ordered response categories.

### 8. Ordered Probit

```r
fit_ordered_probit = econ_ordered(
    Sat ~ Infl + Type + Cont,
    data = housing,
    weights = Freq,
    link = "probit",
    vcov = "HC1"
)

summary(fit_ordered_probit)
```

### 9. Clustered Ordered Logit

Cluster-robust covariance estimation is also supported for ordered models.

```r
housing$cluster_id = rep(
    seq_len(12),
    length.out = nrow(housing)
)

fit_ordered_cluster = econ_ordered(
    Sat ~ Infl + Type + Cont,
    data = housing,
    weights = Freq,
    link = "logit",
    vcov = "CR1",
    cluster = cluster_id
)

summary(fit_ordered_cluster)
```

The cluster identifiers above are illustrative and should not be interpreted as meaningful clusters in the original dataset. In empirical applications, clustering should reflect the sampling design or dependence structure.

### 10. Additional Significance Threshold

```r
# Default significance symbols
summary(fit_ordered_logit)

# Include the additional p < 0.001 threshold
summary(
    fit_ordered_logit,
    extra.significance = TRUE
)
```

### 11. Compatibility with broom and modelsummary

```r
library(broom)
library(modelsummary)

# Coefficient table
tidy(fit_logit)

# Model-level statistics
glance(fit_logit)

# Ordered-model coefficients and thresholds
tidy(fit_ordered_logit)

# Ordered-model confidence intervals
tidy(
    fit_ordered_logit,
    conf.int = TRUE
)

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

The significance symbols used by `modelsummary` are controlled by its own `stars` argument and are independent of the optional `extra.significance` argument in `EconR` summaries.

### 12. Direct Sandwich Compatibility

For model families supporting the corresponding sandwich interfaces, covariance matrices can be computed directly:

```r
sandwich::vcovHC(
    fit_logit,
    type = "HC3"
)
```

## Validation

Version `0.1.3` has been tested against established R implementations and infrastructure, including:

- `lm()`
- `glm()`
- `MASS::polr()`
- `estimatr`
- `sandwich`
- `clubSandwich`
- `broom`

Validation covers:

- Coefficient estimation
- Classical and robust covariance matrices
- Cluster-robust covariance estimation
- Weighted estimation
- Predictions
- Missing-data handling
- Statistical inference
- Model-reporting compatibility

The automated test suite currently contains:

| Model family | Passing tests |
|---|---:|
| Linear models (OLS/WLS) | 490 |
| Binary models (Logit/Probit) | 412 |
| Ordered models (Logit/Probit) | 155 |
| **Total** | **864** |

The complete test suite passes:

```text
FAIL 0 | WARN 0 | SKIP 0 | PASS 1057
```

The package also successfully passes `R CMD check`:

```text
0 errors | 0 warnings | 0 notes
```

These results describe the validation performed for version `0.1.3` and do not imply that all possible model specifications or data configurations have been tested.

## Design Principle

The central design principle of `EconR` is that the estimation sample should be determined by the model being estimated.

An observation is excluded because information required by that model is missing—not because an unrelated variable elsewhere in the data contains a missing value.

`EconR` therefore owns model construction, estimation-sample handling, and estimation semantics, while established packages provide specialized functionality where appropriate.

The implementation combines:

- Consistent model interfaces
- Model-specific missing-data handling
- Explicit covariance-estimator selection
- Standardized inference and reporting
- Compatibility with existing R econometric infrastructure

This philosophy is applied across linear, binary-response, and ordered-response models.

## Development Status

`EconR` is under active development.

Version `0.1.3` extends the validated OLS/WLS and binary Logit/Probit framework with:

- Ordered Logit and Probit estimation
- Robust and clustered covariance estimation for ordered-response models
- Ordered-model coefficient and threshold inference
- Standardized significance symbols across model families
- An optional additional significance threshold
- `broom` integration for ordered models
- An expanded automated test suite with 767 passing tests

Future development will extend the interface to additional econometric model classes.

## Author

David Bustamante Lazo

## Repository

https://github.com/DavidBustamanteL/EconR
