
# EconR

`EconR` is a developing R package for econometric estimation with a consistent interface and model-specific missing-data handling.

The package aims to provide a unified framework for econometric models while preserving familiar R formula syntax and compatibility with established R infrastructure.

## Current Version

**EconR 0.1.3**

The current version implements three model families:

- Ordinary Least Squares (OLS) and Weighted Least Squares (WLS)
- Binary Logit and Probit
- Ordered Logit and Probit

The package supports classical, heteroskedasticity-robust, and cluster-robust covariance estimation, prediction, inference, and integration with standard R model-reporting tools.

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
- Coefficient inference and confidence intervals
- Integration with standard R and `broom` methods

The dependent variable must be numeric and binary, containing only `0` and `1`.

### Ordered Response Models

The `econ_ordered()` function supports:

- Ordered Logit (`link = "logit"`)
- Ordered Probit (`link = "probit"`)
- Weighted maximum-likelihood estimation
- R formula syntax
- Factors, transformations, and interactions
- Subsetting and model-specific missing-data handling
- Classical covariance estimation
- Heteroskedasticity-robust covariance estimation (`HC0`, `HC1`)
- Cluster-robust covariance estimation (`CR0`, `CR1`, `CR1p`, `CR1S`, `stata`)
- Regression coefficient and threshold estimation
- Asymptotic normal inference for classical and HC covariance estimators
- Cluster-adjusted t inference with `G - 1` degrees of freedom
- Wald confidence intervals for coefficients and thresholds
- In-sample and out-of-sample category predictions
- Predicted category probabilities
- Delta-method standard errors for predicted category probabilities
- Weighted and factor-based prediction standard errors
- Response and Pearson residuals
- `na.exclude` restoration of omitted observations
- Log-likelihood, deviance, AIC, and BIC
- Likelihood-ratio tests for nested models using `anova()`
- Integration with `broom::tidy()`, `broom::glance()`, and `broom::augment()`

The dependent variable must be an ordered factor with at least three levels.

Ordered models are estimated using `MASS::polr()`, while `EconR` provides its own covariance estimation, inference, residuals, augmentation, and reporting methods.



#### Ordered-Model Prediction Standard Errors

For an ordered-response model, the probability of category $j$ is

```math
P(Y_i = j \mid X_i) = F(\zeta_j - X_i^\top \beta) - F(\zeta_{j-1} - X_i^\top \beta)
```

where $F$ is the logistic or standard normal cumulative distribution function, $\beta$ contains the regression coefficients, and $\zeta$ contains the category thresholds.

`EconR` uses the **delta method** to estimate standard errors of predicted category probabilities:

```math
SE(\hat{p}_{ij}) = \sqrt{\nabla_\theta p_{ij}^\top \widehat{V} \nabla_\theta p_{ij}}
```

where $\theta = (\beta^\top, \zeta^\top)^\top$ and $\widehat{V}$ is the selected model covariance matrix.

Consequently, prediction standard errors reflect the selected classical, heteroskedasticity-robust, or cluster-robust covariance estimator.

Prediction standard errors have been validated against independent numerical derivatives for both ordered Logit and Probit models, including weighted estimation and factor predictors.


### Missing-Data Handling

`EconR` constructs the estimation sample from the variables actually required by the specified model.

Missing values in unrelated columns therefore do not remove observations from the estimation sample.

Supported `NA` behavior includes:

- `na.omit`
- `na.exclude`
- `na.fail`

Missing values in model inputs such as weights and cluster identifiers are handled within the model-specific estimation workflow.

For ordered-response models, `na.exclude` also preserves exclusion information and restores omitted observations as missing values in supported in-sample predictions and residuals.

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
| `"CR1p"` | Yes | Yes | Yes |
| `"CR1S"` | Yes | Yes | Yes |
| `"CR2"` | Yes | Yes | — |
| `"CR3"` | Yes | Yes | — |
| `"stata"` | Yes | Yes | Yes |

For linear and binary-response models, `EconR` supports all seven heteroskedasticity-consistent covariance estimators available through `sandwich::vcovHC()` and all six cluster-robust covariance estimators available through `clubSandwich::vcovCR()`.

HC estimators account for heteroskedasticity, while CR estimators additionally allow within-cluster error dependence.

The `"stata"` option provides Stata-style cluster-robust covariance estimation.

Ordered-response models implement their supported robust and clustered covariance estimators using analytical observation-level score contributions.

#### Inference Conventions

For binary Logit and Probit models, the reference distribution depends on the covariance estimator:

- **Classical and HC0–HC5:** Asymptotic normal inference (z-tests).
- **CR0, CR1, CR1p, CR1S, CR3, and stata:** Student's t inference with `G - 1` degrees of freedom, where `G` is the number of clusters.
- **CR2:** Student's t inference with coefficient-specific Satterthwaite degrees of freedom, calculated using `clubSandwich`.

For ordered Logit and Probit models:

- **Classical, HC0, and HC1:** Asymptotic normal inference (z-tests).
- **CR0, CR1, CR1p, CR1S, and stata:** Student's t inference with `G - 1` degrees of freedom.

These conventions are applied to coefficient p-values and confidence intervals in `summary()`, `tidy()`, and `confint()`.

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

Ordered-response functionality includes:

- `confint()`
- `anova()`
- `fitted()`
- `residuals()`
- `model.frame()`
- `logLik()`
- `AIC()`
- `BIC()`
- `df.residual()`

Ordered-response models also retain compatibility with relevant methods inherited from `MASS::polr()`.

Compatibility is provided for established reporting and covariance-estimation infrastructure, including:

- `sandwich`
- `broom`
- `modelsummary`

The package supports `broom::tidy()` for coefficient reporting, `broom::glance()` for model-level statistics, and `broom::augment()` for observation-level output.

For ordered models, `augment()` can return category probabilities, predicted categories, ordinal-score residuals for in-sample observations, and delta-method standard errors of category probabilities.

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

# Predicted category probabilities
predict(
    fit_ordered_logit,
    type = "probs"
)

# Predicted categories
predict(
    fit_ordered_logit,
    type = "class"
)
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
    vcov = "CR1S",
    cluster = cluster_id
)

summary(fit_ordered_cluster)

confint(fit_ordered_cluster)
```

The cluster identifiers above are illustrative and should not be interpreted as meaningful clusters in the original dataset. In empirical applications, clustering should reflect the sampling design or dependence structure.

### 10. Ordered Predictions with Standard Errors

`broom::augment()` provides observation-level predictions and delta-method standard errors for category probabilities.

```r
library(broom)

new_housing = housing[1:5, c("Infl", "Type", "Cont")]

ordered_predictions = augment(
    fit_ordered_logit,
    newdata = new_housing,
    se_fit = TRUE
)

ordered_predictions
```

The resulting tibble includes:

- `.fitted`: Most likely response category
- `.fitted_<category>`: Predicted probability of each response category
- `.se.fit_<category>`: Delta-method standard error for each category probability

Category names are converted into valid R column-name components using `make.names()`.

The standard errors use the covariance estimator specified when fitting the model.

### 11. Ordered Residuals

Ordered models provide response and Pearson residuals based on ordinal category scores.

```r
# Ordinal-score response residuals
residuals(
    fit_ordered_logit,
    type = "response"
)

# Standardized ordinal-score residuals
residuals(
    fit_ordered_logit,
    type = "pearson"
)
```

Response residuals are defined as the observed category score minus its predicted expected category score.

Pearson residuals additionally divide by the square root of the predicted variance of the category score.

These are ordinal-score residuals rather than latent-variable residuals.

### 12. Likelihood-Ratio Tests for Ordered Models

Nested ordered-response models can be compared using likelihood-ratio tests.

```r
ordered_small = econ_ordered(
    Sat ~ Infl,
    data = housing,
    weights = Freq,
    link = "logit"
)

ordered_large = econ_ordered(
    Sat ~ Infl + Type + Cont,
    data = housing,
    weights = Freq,
    link = "logit"
)

anova(
    ordered_small,
    ordered_large,
    test = "Chisq"
)
```

Models must use the same link function, estimation sample, dependent variable, and observation weights. Their specifications must be nested and supplied in increasing order of complexity.

### 13. Additional Significance Threshold

```r
# Default significance symbols
summary(fit_ordered_logit)

# Include the additional p < 0.001 threshold
summary(
    fit_ordered_logit,
    extra.significance = TRUE
)
```

### 14. Compatibility with broom and modelsummary

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

# Ordered-model fit statistics
glance(fit_ordered_logit)

# Ordered-model observation-level predictions
augment(
    fit_ordered_logit,
    se_fit = TRUE
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

### 15. Direct Sandwich Compatibility

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
- Delta-method prediction standard errors
- Factor-based prediction
- Missing-data handling and `na.exclude` restoration
- Statistical inference
- Model-reporting compatibility

The automated test suite currently reports:

| Model family | Passing assertions |
|---|---:|
| Linear models (OLS/WLS) | 490 |
| Binary models (Logit/Probit) | 412 |
| Ordered models (Logit/Probit) | 443 |
| **Total** | **1,345** |

The complete test suite passes:

```text
FAIL 0 | WARN 0 | SKIP 0 | PASS 1345
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
- Delta-method prediction standard errors
- Weighted and factor-based ordered predictions
- Ordered-model residuals and `na.exclude` restoration
- Likelihood-ratio comparisons of nested ordered models
- Standardized significance symbols across model families
- An optional additional significance threshold
- Expanded `broom` integration for ordered models
- An automated test suite with 1,345 passing assertions
- A clean `R CMD check` with no errors, warnings, or notes

Future development will extend the interface to additional econometric model classes.

## Author

David Bustamante Lazo

## Repository

https://github.com/DavidBustamanteL/EconR
