
# EconR

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23547786.svg)](https://doi.org/10.5281/zenodo.23547786)

`EconR` is an R package for econometric estimation with a consistent model interface, model-specific missing-data handling, and flexible covariance estimation.

It combines familiar R formula syntax with standardized inference, prediction, and reporting across linear, binary, ordered, and multinomial regression models.

**Current version: 0.1.5**

## 1. Supported Models

| Function | Model | Estimation |
|---|---|---|
| `econ_lm()` | Linear regression | OLS and WLS |
| `econ_binary()` | Binary Logit and Probit | Maximum likelihood |
| `econ_ordered()` | Ordered Logit and Probit | Maximum likelihood |
| `econ_multinomial()` | Multinomial Logit and Probit | Maximum likelihood |

All models support observation weights, subsetting, and model-specific missing-data handling.

### Linear Regression

`econ_lm()` supports OLS and WLS with:

- Factors, interactions, transformations, and custom contrasts
- Offsets and rank-deficient specifications
- Classical, heteroskedasticity-robust, and cluster-robust covariance estimation
- Standard inference, prediction, and model diagnostics

### Binary Logit and Probit

`econ_binary()` estimates binary-response models using:

- `link = "logit"` or `link = "probit"`
- A numeric response coded `0` and `1`
- Classical, HC, and cluster-robust covariance estimators
- Predicted probabilities, linear predictors, and confidence intervals
- Likelihood-based goodness-of-fit statistics

### Ordered Logit and Probit

`econ_ordered()` estimates models with an **ordered factor response containing at least three levels**.

Features include:

- `link = "logit"` or `link = "probit"`
- Regression coefficients and category thresholds
- Classical, HC0/HC1, and supported cluster-robust covariance estimators
- Predicted categories and category probabilities
- Delta-method standard errors for predicted probabilities
- Response and Pearson residuals based on ordinal scores
- Likelihood-ratio comparisons of nested models

Estimation uses `MASS::polr()`, with additional covariance, inference, prediction, and reporting functionality implemented by EconR.

### Multinomial Logit and Probit

`econ_multinomial()` estimates models for an **unordered factor response containing at least three observed categories**.

The first response level is the reference category.

| Specification | Arguments |
|---|---|
| Multinomial Logit | `link = "logit"` |
| Multinomial Probit, fixed error covariance | `link = "probit", error_covariance = "fixed"` |
| Multinomial Probit, estimated error covariance | `link = "probit", error_covariance = "estimated"` |

Multinomial Logit uses `nnet::multinom()`. Multinomial Probit uses EconR's multivariate-normal probability and maximum-likelihood implementation.

Supported functionality includes:

- Classical, HC0/HC1, and cluster-robust covariance estimation
- Coefficient inference and joint Wald tests
- Predicted categories and choice probabilities
- Delta-method prediction standard errors
- Average marginal effects for numeric predictors
- Average discrete probability changes for factor predictors
- Log-likelihood, AIC, BIC, and McFadden's pseudo-R-squared
- `broom` integration for model and observation-level reporting

For Multinomial Probit, estimating the error covariance introduces additional parameters and increases computational cost.

## 2. Covariance Estimation

The available covariance estimators depend on the model family.

| Estimator | Linear | Binary | Ordered | Multinomial |
|---|:---:|:---:|:---:|:---:|
| `classical` | Yes | Yes | Yes | Yes |
| `HC0`, `HC1` | Yes | Yes | Yes | Yes |
| `HC2`, `HC3` | Yes | Yes | — | — |
| `HC4`, `HC4m`, `HC5` | Yes | Yes | — | — |
| `CR0`, `CR1` | Yes | Yes | Yes | Yes |
| `CR1p`, `CR1S` | Yes | Yes | Yes | Yes |
| `CR2`, `CR3` | Yes | Yes | — | — |
| `stata` | Yes | Yes | Yes | Yes |

Select the covariance estimator directly when fitting a model:

```r
fit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "HC3"
)

summary(fit)
vcov(fit)
```

For cluster-robust estimation, specify both `vcov` and `cluster`.

Classical and HC inference uses asymptotic normal distributions for nonlinear models. Supported cluster-robust inference uses cluster-adjusted t distributions, with Satterthwaite degrees of freedom for CR2 in binary models.

For detailed covariance definitions, see the relevant function help pages.

## 3. Missing-Data Handling

EconR constructs the estimation sample from the variables required by the specified model.

Missing values in unrelated dataset columns do not automatically exclude observations.

Supported missing-data actions include:

- `na.omit`
- `na.exclude`
- `na.fail`

Model inputs such as observation weights and cluster identifiers are incorporated into the relevant estimation workflow.

Where supported, `na.exclude` preserves the original observation alignment when returning in-sample predictions and residuals.

## 4. Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("DavidBustamanteL/EconR")
```

To install the tagged 0.1.5 release, once the tag is published:

```r
remotes::install_github(
    "DavidBustamanteL/EconR@v0.1.5"
)
```

Load the package:

```r
library(EconR)
```

## 5. Quick Examples

### 5.1. Linear Regression

```r
fit_lm = econ_lm(
    mpg ~ wt * factor(am) + hp,
    data = mtcars,
    vcov = "HC3"
)

summary(fit_lm)
confint(fit_lm)
predict(fit_lm, newdata = mtcars[1:3, ])
```

For cluster-robust estimation:

```r
fit_cluster = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    vcov = "CR2",
    cluster = cyl
)

summary(fit_cluster)
```

The cluster variable in this example is illustrative; empirical clustering should follow the data's dependence or sampling structure.

### 5.2. Binary Logit and Probit

```r
fit_logit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "HC3"
)

fit_probit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "probit",
    vcov = "HC1"
)

summary(fit_logit)
summary(fit_probit)

# Predicted probabilities
predict(fit_logit, type = "response")
```

### 5.3. Ordered Logit and Probit

The `MASS::housing` dataset contains an ordered response (`Sat`) and observation frequencies (`Freq`).

```r
library(MASS)

fit_ordered_logit = econ_ordered(
    Sat ~ Infl + Type + Cont,
    data = housing,
    weights = Freq,
    link = "logit",
    vcov = "HC1"
)

fit_ordered_probit = econ_ordered(
    Sat ~ Infl + Type + Cont,
    data = housing,
    weights = Freq,
    link = "probit"
)

summary(fit_ordered_logit)
summary(fit_ordered_probit)

# Category probabilities
predict(fit_ordered_logit, type = "probs")

# Most likely category
predict(fit_ordered_logit, type = "class")
```

For prediction standard errors:

```r
broom::augment(
    fit_ordered_logit,
    newdata = housing[1:5, c("Infl", "Type", "Cont")],
    se_fit = TRUE
)
```

This returns category predictions, probabilities, and delta-method standard errors.

### 5.4. Multinomial Logit and Probit

The following example creates a three-category response with `A` as the reference category.

```r
set.seed(123)

n = 200

dat = data.frame(
    x = rnorm(n)
)

eta_B = 0.5 + 0.8 * dat$x
eta_C = -0.3 - 0.5 * dat$x

P = cbind(
    A = 1,
    B = exp(eta_B),
    C = exp(eta_C)
)

P = P / rowSums(P)

dat$y = factor(
    vapply(
        seq_len(n),
        function(i) {
            sample(
                c("A", "B", "C"),
                size = 1,
                prob = P[i, ]
            )
        },
        character(1)
    ),
    levels = c("A", "B", "C")
)
```

**Multinomial Logit:**

```r
fit_multi_logit = econ_multinomial(
    y ~ x,
    data = dat,
    link = "logit"
)

summary(fit_multi_logit)

predict(
    fit_multi_logit,
    newdata = dat[1:5, ],
    type = "probs"
)
```

**Multinomial Probit with fixed error covariance:**

```r
fit_multi_probit = econ_multinomial(
    y ~ x,
    data = dat,
    link = "probit",
    error_covariance = "fixed"
)

summary(fit_multi_probit)

predict(
    fit_multi_probit,
    newdata = dat[1:5, ],
    type = "probs"
)
```

To estimate the Probit error covariance instead, set:

```r
error_covariance = "estimated"
```

inside the `econ_multinomial()` call.

### 5.5. Multinomial Marginal Effects

Marginal effects express changes in predicted choice probabilities rather than changes in log-odds or latent utilities.

```r
marginal_effects(
    fit_multi_logit
)
```

For numeric predictors, EconR computes average marginal effects. For factor predictors, it computes average discrete changes in category probabilities.

Standard errors use the delta method and the covariance estimator selected during estimation.

### 5.6. Model Reporting

EconR supports standard model methods, including:

- `summary()`, `coef()`, and `vcov()`
- `predict()`, `confint()`, and `logLik()` where implemented
- `AIC()` and `BIC()` for likelihood-based models
- `broom::tidy()`, `broom::glance()`, and `broom::augment()`

For example:

```r
broom::tidy(fit_multi_logit)
broom::glance(fit_multi_logit)

broom::augment(
    fit_multi_logit,
    newdata = dat[1:5, ],
    se_fit = TRUE
)
```

Regression tables can also be generated using `modelsummary`:

```r
modelsummary::modelsummary(
    list(
        "Binary Logit" = fit_logit,
        "Binary Probit" = fit_probit
    )
)
```

## 6. Validation

EconR has been validated against established R implementations and econometric infrastructure, including:

- `lm()` and `estimatr` for linear regression
- `glm()` for binary-response models
- `MASS::polr()` for ordered-response models
- `nnet` and `mlogit` for multinomial models
- `sandwich` and `clubSandwich` for covariance estimation
- `broom` for reporting compatibility

Validation covers parameter estimation, covariance matrices, inference, observation weights, missing-data handling, predictions, and relevant reporting methods.

The latest reported model-specific test results are:

| Model family | Passing test expectations |
|---|---:|
| Linear (OLS/WLS) | 490 |
| Binary Logit/Probit | 424 |
| Ordered Logit/Probit | 467 |
| Multinomial Logit/Probit | 485 |
| **Total** | **1,866** |

The multinomial tests include fixed- and estimated-covariance Probit specifications.

An additional console comparison of Multinomial Logit against `mlogit` confirmed close agreement in coefficients, standard errors, log-likelihood, AIC, and predicted probabilities.

The complete package also passed `R CMD check` on Windows with R 4.6.1:

```text
0 errors | 0 warnings | 0 notes
```

These results document the tested specifications and do not imply that every possible model configuration has been validated.

## 7. Design Principle

EconR's central principle is that **the model determines the estimation sample**.

An observation should be excluded because information required by the model is missing, not because an unrelated variable in the dataset contains an `NA`.

The package combines this principle with:

- Consistent function interfaces
- Explicit covariance-estimator selection
- Model-appropriate inference
- Standard R prediction and reporting methods

EconR is under active development.

Version 0.1.5 introduces Multinomial Logit and Probit, including fixed and estimated Probit error covariance specifications, robust and clustered inference, category predictions, and marginal effects.

## 8. Documentation

For function-specific arguments and supported methods:

```r
help(EconR)

?econ_lm
?econ_binary
?econ_ordered
?econ_multinomial
```

## 9. Author and Citation

**Author:** David Bustamante Lazo

**Repository:** https://github.com/DavidBustamanteL/EconR

**DOI:** https://doi.org/10.5281/zenodo.23547786

If you use EconR in academic research, please cite:

Bustamante Lazo, D. (2026). *EconR: Econometric Models with Consistent Estimation Semantics* (Version 0.1.5) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23547786
