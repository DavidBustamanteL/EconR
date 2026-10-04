################################ Testing Binary Models EconR ###############################
############################################################################################

# let us begin as usual by setting our wd
rm(list = ls())
gc()

# English Format and Avoiding scientific notation
Sys.setlocale("LC_TIME", "English")
options(scipen = 999)

# Packages
library(tidyverse)
library(magrittr)
library(estimatr)
library(sandwich)
library(clubSandwich)

# Installing test pack if required
if (!requireNamespace("EconR", quietly = TRUE)) {
    install.packages(
        "G:/My Drive/5. CURSOS ONLINE/R Practice FULL/Econ Own Package",
        repos = NULL,
        type = "source"
    )
}

library(EconR)


#### 1. Basic Logit ####

## With mtcars ##
econ_logit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit"
)

glm_logit = glm(
    am ~ wt + hp,
    data = mtcars,
    family = binomial(
        link = "logit"
    )
)

coef(econ_logit)

coef(glm_logit)

all.equal(
    coef(econ_logit),
    coef(glm_logit)
)


#### 2. Logit Fitted Probabilities ####

econ_fitted = econ_logit$fitted.values

glm_fitted = fitted(
    glm_logit
)

cbind(
    econ_fitted,
    glm_fitted,
    difference = econ_fitted - glm_fitted
)

all.equal(
    econ_fitted,
    glm_fitted
)


#### 3. Logit Linear Predictors ####

econ_link = econ_logit$linear.predictors

glm_link = predict(
    glm_logit,
    type = "link"
)

all.equal(
    econ_link,
    glm_link
)


#### 4. Logit Deviance ####

comparison = data.frame(
    statistic = c(
        "Residual deviance",
        "Null deviance"
    ),
    EconR = c(
        econ_logit$deviance,
        econ_logit$null.deviance
    ),
    GLM = c(
        glm_logit$deviance,
        glm_logit$null.deviance
    )
)

print(comparison)

all.equal(
    econ_logit$deviance,
    glm_logit$deviance
)

all.equal(
    econ_logit$null.deviance,
    glm_logit$null.deviance
)


#### 5. Classical Logit Covariance ####

glm_vcov = vcov(
    glm_logit
)

X = econ_logit$x

w = econ_logit$weights

p = econ_logit$fitted.values

W = w * p * (1 - p)

econ_vcov_manual = solve(
    crossprod(
        X,
        X * W
    )
)

print(glm_vcov)

print(econ_vcov_manual)

all.equal(
    unname(econ_vcov_manual),
    unname(glm_vcov),
    tolerance = 1e-8
)


#### 5.1. Classical Covariance via QR ####

qr_fit = glm_logit$qr

p = glm_logit$rank

R = qr.R(
    qr_fit
)

R = R[
    seq_len(p),
    seq_len(p),
    drop = FALSE
]

V_qr = chol2inv(
    R
)

dimnames(V_qr) = dimnames(
    glm_vcov
)

print(V_qr)

all.equal(
    V_qr,
    glm_vcov,
    tolerance = 1e-10
)


#### 6. EconR Classical Covariance Method ####

econ_logit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit"
)

all.equal(
    vcov(econ_logit),
    vcov(glm_logit),
    tolerance = 1e-10
)


#### 7. Classical Logit Inference ####

V = vcov(
    econ_logit
)

estimate = coef(
    econ_logit
)

std_error = sqrt(
    diag(V)
)

z_value = estimate / std_error

p_value = 2 * stats::pnorm(
    abs(z_value),
    lower.tail = FALSE
)

econ_table = cbind(
    Estimate = estimate,
    "Std. Error" = std_error,
    "z value" = z_value,
    "Pr(>|z|)" = p_value
)

glm_table = summary(
    glm_logit
)$coefficients

print(econ_table)

all.equal(
    econ_table,
    glm_table,
    tolerance = 1e-10
)


#### 8. Logit Summary ####

econ_logit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit"
)

summary(econ_logit)

all.equal(
    summary(econ_logit)$coefficients,
    summary(glm_logit)$coefficients,
    tolerance = 1e-10
)


#### 9. Logit Information Criteria ####

eta = econ_logit$linear.predictors

y = econ_logit$y

w = econ_logit$weights

log_p = -log1p(
    exp(-eta)
)

log_1mp = -log1p(
    exp(eta)
)

loglik_econ = sum(
    w * (
        y * log_p +
            (1 - y) * log_1mp
    )
)

k = econ_logit$rank

n = sum(w > 0)

aic_econ = -2 * loglik_econ + 2 * k

bic_econ = -2 * loglik_econ + log(n) * k

comparison = cbind(
    EconR = c(
        loglik_econ,
        aic_econ,
        bic_econ
    ),
    GLM = c(
        as.numeric(logLik(glm_logit)),
        AIC(glm_logit),
        BIC(glm_logit)
    )
)

rownames(comparison) = c(
    "LogLik",
    "AIC",
    "BIC"
)

print(comparison)

all.equal(
    unname(comparison[, "EconR"]),
    unname(comparison[, "GLM"]),
    tolerance = 1e-10
)


#### 10. Log-Likelihood S3 Compatibility ####

devtools::load_all()

econ_logit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit"
)

all.equal(
    as.numeric(logLik(econ_logit)),
    as.numeric(logLik(glm_logit)),
    tolerance = 1e-10
)

all.equal(
    AIC(econ_logit),
    AIC(glm_logit),
    tolerance = 1e-10
)

all.equal(
    BIC(econ_logit),
    BIC(glm_logit),
    tolerance = 1e-10
)


#### 11. Logit In-Sample Prediction ####

econ_logit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit"
)

all.equal(
    predict(econ_logit, type = "link"),
    predict(glm_logit, type = "link")
)

all.equal(
    predict(econ_logit, type = "response"),
    predict(glm_logit, type = "response")
)

#### 12. Logit Out-of-Sample Prediction ####

new_data = data.frame(
    wt = c(2.0, 2.5, 3.0, 3.5, 4.0),
    hp = c(100, 120, 150, 180, 200)
)

# 12.1. Linear predictors #
econ_link = predict(
    econ_logit,
    newdata = new_data,
    type = "link"
)

glm_link = predict(
    glm_logit,
    newdata = new_data,
    type = "link"
)

all.equal(
    econ_link,
    glm_link
)

# 12.2. Predicted probabilities #
econ_prob = predict(
    econ_logit,
    newdata = new_data,
    type = "response"
)

glm_prob = predict(
    glm_logit,
    newdata = new_data,
    type = "response"
)

all.equal(
    econ_prob,
    glm_prob
)


#### 13. Logit Prediction with Factors ####

# 13.1. Prepare data #
factor_data = mtcars
factor_data$cyl = factor(factor_data$cyl)

# 13.2. Estimate models #
econ_factor = econ_binary(
    am ~ wt + cyl,
    data = factor_data,
    link = "logit"
)

glm_factor = glm(
    am ~ wt + cyl,
    data = factor_data,
    family = binomial(link = "logit")
)

# 13.3. Construct new observations #
new_factor_data = data.frame(
    wt = c(2.5, 3.0, 3.5),
    cyl = factor(
        c(4, 6, 8),
        levels = levels(factor_data$cyl)
    )
)

# 13.4. Compare predictions #
all.equal(
    predict(econ_factor, newdata = new_factor_data, type = "link"),
    predict(glm_factor, newdata = new_factor_data, type = "link")
)

all.equal(
    predict(econ_factor, newdata = new_factor_data, type = "response"),
    predict(glm_factor, newdata = new_factor_data, type = "response")
)


#### 14. Prediction with Missing Values ####

# 14.1. Construct new observations #
new_na_data = data.frame(
    wt = c(2.0, NA, 3.0, 3.5),
    hp = c(100, 120, NA, 180)
)

# 14.2. Compare linear predictions #
all.equal(
    predict(
        econ_logit,
        newdata = new_na_data,
        type = "link"
    ),
    predict(
        glm_logit,
        newdata = new_na_data,
        type = "link"
    )
)

# 14.3. Compare predicted probabilities #
all.equal(
    predict(
        econ_logit,
        newdata = new_na_data,
        type = "response"
    ),
    predict(
        glm_logit,
        newdata = new_na_data,
        type = "response"
    )
)


#### 15. Weighted Logit #### -> warning expected

# 15.1. Prepare weights #
weighted_data = mtcars

weighted_data$obs_weight = seq(
    0.5,
    2,
    length.out = nrow(weighted_data)
)

# 15.2. Estimate models #
econ_weighted = econ_binary(
    am ~ wt + hp,
    data = weighted_data,
    weights = obs_weight,
    link = "logit"
)

glm_weighted = glm(
    am ~ wt + hp,
    data = weighted_data,
    weights = obs_weight,
    family = binomial(link = "logit")
)

# 15.3. Compare coefficients #
all.equal(
    coef(econ_weighted),
    coef(glm_weighted)
)

# 15.4. Compare covariance matrices #
all.equal(
    vcov(econ_weighted),
    vcov(glm_weighted),
    tolerance = 1e-10
)

# 15.5. Compare fitted probabilities #
all.equal(
    econ_weighted$fitted.values,
    fitted(glm_weighted)
)


#### 16. Logit Subset Estimation ####

# 16.1. Estimate models #
econ_subset = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    subset = cyl != 6,
    link = "logit"
)

glm_subset = glm(
    am ~ wt + hp,
    data = mtcars,
    subset = cyl != 6,
    family = binomial(link = "logit")
)

# 16.2. Compare coefficients #
all.equal(
    coef(econ_subset),
    coef(glm_subset)
)

# 16.3. Compare covariance matrices #
all.equal(
    vcov(econ_subset),
    vcov(glm_subset),
    tolerance = 1e-10
)

# 16.4. Compare estimation sample #
all.equal(
    nrow(econ_subset$model),
    nobs(glm_subset)
)


#### 17. Logit Estimation with Missing Values ####

# 17.1. Prepare data #
na_data = mtcars

na_data$wt[c(2, 5)] = NA
na_data$hp[c(7, 12)] = NA

# 17.2. Estimate models #
econ_na = econ_binary(
    am ~ wt + hp,
    data = na_data,
    link = "logit"
)

glm_na = glm(
    am ~ wt + hp,
    data = na_data,
    family = binomial(link = "logit")
)

# 17.3. Compare coefficients #
all.equal(
    coef(econ_na),
    coef(glm_na)
)

# 17.4. Compare covariance matrices #
all.equal(
    vcov(econ_na),
    vcov(glm_na),
    tolerance = 1e-10
)

# 17.5. Compare estimation sample #
all.equal(
    nrow(econ_na$model),
    nobs(glm_na)
)

# 17.6. Check omitted observations #
econ_na$omitted


#### 18. Logit with Missing Observation Weights ####

# 18.1. Prepare data #
weight_na_data = mtcars

weight_na_data$obs_weight = rep(
    1,
    nrow(weight_na_data)
)

weight_na_data$obs_weight[c(3, 9)] = NA

# 18.2. Estimate models #
econ_weight_na = econ_binary(
    am ~ wt + hp,
    data = weight_na_data,
    weights = obs_weight,
    link = "logit"
)

glm_weight_na = glm(
    am ~ wt + hp,
    data = weight_na_data,
    weights = obs_weight,
    family = binomial(link = "logit")
)

# 18.3. Compare coefficients #
all.equal(
    coef(econ_weight_na),
    coef(glm_weight_na)
)

# 18.4. Compare covariance matrices #
all.equal(
    vcov(econ_weight_na),
    vcov(glm_weight_na),
    tolerance = 1e-10
)

# 18.5. Compare sample sizes #
all.equal(
    nrow(econ_weight_na$model),
    nobs(glm_weight_na)
)

# 18.6. Inspect omitted observations #
econ_weight_na$omitted


#### 19. Probit Estimation ####

# 19.1. Estimate models #
econ_probit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "probit"
)

glm_probit = glm(
    am ~ wt + hp,
    data = mtcars,
    family = binomial(link = "probit")
)

# 19.2. Compare coefficients #
all.equal(
    coef(econ_probit),
    coef(glm_probit)
)

# 19.3. Compare covariance matrices #
all.equal(
    vcov(econ_probit),
    vcov(glm_probit),
    tolerance = 1e-10
)

# 19.4. Compare fitted probabilities #
all.equal(
    econ_probit$fitted.values,
    fitted(glm_probit)
)


#### 20. Probit Prediction and Inference ####

# 20.1. Compare coefficient inference #
all.equal(
    summary(econ_probit)$coefficients,
    summary(glm_probit)$coefficients,
    tolerance = 1e-10
)

# 20.2. Compare in-sample probabilities #
all.equal(
    predict(econ_probit, type = "response"),
    predict(glm_probit, type = "response")
)

# 20.3. Compare out-of-sample linear predictions #
all.equal(
    predict(
        econ_probit,
        newdata = new_data,
        type = "link"
    ),
    predict(
        glm_probit,
        newdata = new_data,
        type = "link"
    )
)

# 20.4. Compare out-of-sample probabilities #
all.equal(
    predict(
        econ_probit,
        newdata = new_data,
        type = "response"
    ),
    predict(
        glm_probit,
        newdata = new_data,
        type = "response"
    )
)


#### 21. Probit Information Criteria ####

# 21.1. Compare log-likelihood #
all.equal(
    as.numeric(logLik(econ_probit)),
    as.numeric(logLik(glm_probit)),
    tolerance = 1e-10
)

# 21.2. Compare AIC #
all.equal(
    AIC(econ_probit),
    AIC(glm_probit),
    tolerance = 1e-10
)

# 21.3. Compare BIC #
all.equal(
    BIC(econ_probit),
    BIC(glm_probit),
    tolerance = 1e-10
)


#### 22. Weighted Probit ####

# 22.1. Estimate models #
econ_probit_w = econ_binary(
    am ~ wt + hp,
    data = weighted_data,
    weights = obs_weight,
    link = "probit"
)

glm_probit_w = glm(
    am ~ wt + hp,
    data = weighted_data,
    weights = obs_weight,
    family = binomial(link = "probit")
)

# 22.2. Compare coefficients #
all.equal(
    coef(econ_probit_w),
    coef(glm_probit_w)
)

# 22.3. Compare covariance matrices #
all.equal(
    vcov(econ_probit_w),
    vcov(glm_probit_w),
    tolerance = 1e-10
)

# 22.4. Compare fitted probabilities #
all.equal(
    econ_probit_w$fitted.values,
    fitted(glm_probit_w)
)


#### 23. Binary Model Score Contributions ####

# 23.1. Re-estimate logit #
econ_logit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit"
)

# 23.2. Re-estimate probit #
econ_probit = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "probit"
)

# 23.3. Compare score contributions #
all.equal(
    sandwich::estfun(econ_logit),
    sandwich::estfun(glm_logit),
    tolerance = 1e-10
)

all.equal(
    sandwich::estfun(econ_probit),
    sandwich::estfun(glm_probit),
    tolerance = 1e-10
)


#### 24. Binary Model Bread Matrix ####

# 24.1. Logit #
all.equal(
    sandwich::bread(econ_logit),
    sandwich::bread(glm_logit),
    tolerance = 1e-10
)

# 24.2. Probit #
all.equal(
    sandwich::bread(econ_probit),
    sandwich::bread(glm_probit),
    tolerance = 1e-10
)


#### 25. Binary Robust Covariance ####

# 25.1. Logit HC0 #
all.equal(
    sandwich::vcovHC(econ_logit, type = "HC0"),
    sandwich::vcovHC(glm_logit, type = "HC0"),
    tolerance = 1e-10
)

# 25.2. Logit HC1 #
all.equal(
    sandwich::vcovHC(econ_logit, type = "HC1"),
    sandwich::vcovHC(glm_logit, type = "HC1"),
    tolerance = 1e-10
)

# 25.3. Probit HC0 #
all.equal(
    sandwich::vcovHC(econ_probit, type = "HC0"),
    sandwich::vcovHC(glm_probit, type = "HC0"),
    tolerance = 1e-10
)

# 25.4. Probit HC1 #
all.equal(
    sandwich::vcovHC(econ_probit, type = "HC1"),
    sandwich::vcovHC(glm_probit, type = "HC1"),
    tolerance = 1e-10
)


#### 26. Binary Model Leverage ####

# 26.1. Logit #
all.equal(
    hatvalues(econ_logit),
    hatvalues(glm_logit),
    tolerance = 1e-10
)

# 26.2. Probit #
all.equal(
    hatvalues(econ_probit),
    hatvalues(glm_probit),
    tolerance = 1e-10
)


#### 27. Binary HC2 and HC3 Covariance ####

# 27.1. Logit HC2 #
all.equal(
    sandwich::vcovHC(econ_logit, type = "HC2"),
    sandwich::vcovHC(glm_logit, type = "HC2"),
    tolerance = 1e-10
)

# 27.2. Logit HC3 #
all.equal(
    sandwich::vcovHC(econ_logit, type = "HC3"),
    sandwich::vcovHC(glm_logit, type = "HC3"),
    tolerance = 1e-10
)

# 27.3. Probit HC2 #
all.equal(
    sandwich::vcovHC(econ_probit, type = "HC2"),
    sandwich::vcovHC(glm_probit, type = "HC2"),
    tolerance = 1e-10
)

# 27.4. Probit HC3 #
all.equal(
    sandwich::vcovHC(econ_probit, type = "HC3"),
    sandwich::vcovHC(glm_probit, type = "HC3"),
    tolerance = 1e-10
)


#### 28. EconR Robust Covariance Interface ####

# 28.1. Logit #
sapply(
    c("HC0", "HC1", "HC2", "HC3"),
    function(type) {
        isTRUE(all.equal(
            vcov(econ_logit, type = type),
            sandwich::vcovHC(glm_logit, type = type),
            tolerance = 1e-10
        ))
    }
)

# 28.2. Probit #
sapply(
    c("HC0", "HC1", "HC2", "HC3"),
    function(type) {
        isTRUE(all.equal(
            vcov(econ_probit, type = type),
            sandwich::vcovHC(glm_probit, type = type),
            tolerance = 1e-10
        ))
    }
)


#### 29. Persistent Robust Covariance ####

econ_logit_hc3 = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "HC3"
)

# 29.1. Verify stored selection #
econ_logit_hc3$vcov_type

# 29.2. Verify default covariance #
all.equal(
    vcov(econ_logit_hc3),
    sandwich::vcovHC(glm_logit, type = "HC3"),
    tolerance = 1e-10
)

# 29.3. Verify classical override #
all.equal(
    vcov(econ_logit_hc3, type = "classical"),
    vcov(glm_logit),
    tolerance = 1e-10
)


#### 30. Robust Binary Model Summary ####

# 30.1. Extract EconR inference #
econ_summary = summary(econ_logit_hc3)$coefficients

# 30.2. Construct reference inference #
V = sandwich::vcovHC(
    glm_logit,
    type = "HC3"
)

estimate = coef(glm_logit)
std_error = sqrt(diag(V))
z_value = estimate / std_error
p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE)

glm_summary = cbind(
    Estimate = estimate,
    "Std. Error" = std_error,
    "z value" = z_value,
    "Pr(>|z|)" = p_value
)

# 30.3. Compare complete inference #
all.equal(
    econ_summary,
    glm_summary,
    tolerance = 1e-10
)


#### 31. Robust Probit Summary ####

# 31.1. Estimate robust probit #
econ_probit_hc1 = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "probit",
    vcov = "HC1"
)

# 31.2. Construct reference inference #
V = sandwich::vcovHC(
    glm_probit,
    type = "HC1"
)

estimate = coef(glm_probit)
std_error = sqrt(diag(V))
z_value = estimate / std_error
p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE)

glm_summary = cbind(
    Estimate = estimate,
    "Std. Error" = std_error,
    "z value" = z_value,
    "Pr(>|z|)" = p_value
)

# 31.3. Compare complete inference #
all.equal(
    summary(econ_probit_hc1)$coefficients,
    glm_summary,
    tolerance = 1e-10
)


#### 33. Binary Cluster Alignment ####

# 33.1. Estimate clustered model #
econ_cluster = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "CR0",
    cluster = cyl
)

# 33.2. Inspect cluster IDs #
econ_cluster$cluster

# 33.3. Compare with original data #
all.equal(
    as.vector(econ_cluster$cluster),
    mtcars$cyl,
    check.attributes = FALSE
)

# 33.4. Check sample alignment #
length(econ_cluster$cluster) == nrow(econ_cluster$model)


#### 34. Binary CR0 Covariance ####

# 34.1. Estimate reference model #
glm_cluster = glm(
    am ~ wt + hp,
    data = mtcars,
    family = binomial(link = "logit")
)

# 34.2. Construct reference covariance #
V_reference = sandwich::vcovCL(
    glm_cluster,
    cluster = mtcars$cyl,
    type = "HC",
    cadjust = FALSE
)

# 34.3. Compare covariance matrices #
all.equal(
    vcov(econ_cluster),
    V_reference,
    tolerance = 1e-10
)


#### 35. CR0 Cross-Package Validation #### -> tolerance to 1e-6 here

# 35.1. Construct clubSandwich reference #
V_club = clubSandwich::vcovCR(
    glm_cluster,
    cluster = mtcars$cyl,
    type = "CR0"
)

# 35.2. Compare covariance matrices #
all.equal(
    unclass(vcov(econ_cluster)),
    unclass(as.matrix(V_club)),
    tolerance = 1e-6
)


#### 36. Binary CR2 Reference ####

# 36.1. Construct CR2 covariance #
V_CR2_reference = clubSandwich::vcovCR(
    glm_cluster,
    cluster = mtcars$cyl,
    type = "CR2"
)

# 36.2. Inspect result #
V_CR2_reference

# 36.2. Reconstruct GLM from EconR Sample
glm_reconstructed = glm(
    formula = econ_cluster$formula,
    data = econ_cluster$model,
    weights = econ_cluster$weights,
    family = econ_cluster$family
)

# Compare coefficients
all.equal(
    coef(glm_reconstructed),
    coef(glm_cluster),
    tolerance = 1e-10
)

# Compare CR2 covariance
all.equal(
    unclass(as.matrix(clubSandwich::vcovCR(
        glm_reconstructed,
        cluster = econ_cluster$cluster,
        type = "CR2"
    ))),
    unclass(as.matrix(V_CR2_reference)),
    tolerance = 1e-10
)


#### 37. Binary CR2 Covariance ####

econ_CR2 = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "CR2",
    cluster = cyl
)

all.equal(
    unclass(vcov(econ_CR2)),
    unclass(as.matrix(V_CR2_reference)),
    tolerance = 1e-10
)


#### 38. Binary CR2 Sample Alignment ####

# 38.1. Create test data #
df_CR2 = mtcars
df_CR2$hp[c(2, 7)] = NA

# 38.2. Estimate EconR model #
econ_CR2_sample = econ_binary(
    am ~ wt + hp,
    data = df_CR2,
    subset = mpg > 18,
    link = "logit",
    vcov = "CR2",
    cluster = cyl
)

# 38.3. Construct reference using EconR's estimation sample #
glm_CR2_sample = glm(
    am ~ wt + hp,
    data = econ_CR2_sample$model,
    family = binomial(link = "logit")
)

V_CR2_sample_reference = clubSandwich::vcovCR(
    glm_CR2_sample,
    cluster = econ_CR2_sample$cluster,
    type = "CR2"
)

# 38.4. Compare covariance matrices #
all.equal(
    unclass(vcov(econ_CR2_sample)),
    unclass(as.matrix(V_CR2_sample_reference)),
    tolerance = 1e-10
)


#### 39. Binary Stata-Style Clustered Covariance ####

# 39.1. Estimate EconR model #
econ_stata = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "stata",
    cluster = cyl
)

# 39.2. Construct reference covariance #
V_stata_reference = sandwich::vcovCL(
    glm_cluster,
    cluster = mtcars$cyl,
    type = "HC1"
)

# 39.3. Compare covariance matrices #
all.equal(
    vcov(econ_stata),
    V_stata_reference,
    tolerance = 1e-10
)


#### 40. Binary Clustered Inference ####

# 40.1. CR0 standard errors #
all.equal(
    summary(econ_cluster)$coefficients[, 2],
    sqrt(diag(vcov(econ_cluster))),
    tolerance = 1e-10
)

# 40.2. CR2 standard errors #
all.equal(
    summary(econ_CR2)$coefficients[, 2],
    sqrt(diag(vcov(econ_CR2))),
    tolerance = 1e-10
)

# 40.3. Stata-style standard errors #
all.equal(
    summary(econ_stata)$coefficients[, 2],
    sqrt(diag(vcov(econ_stata))),
    tolerance = 1e-10
)


#### 41. broom and modelsummary Compatibility ####

# 41.1. Estimate binary model #
econ_compat = econ_binary(
    am ~ wt + hp,
    data = mtcars,
    link = "logit",
    vcov = "HC3"
)

# 41.2. broom methods #
broom::tidy(econ_compat)
broom::glance(econ_compat)

## 41.3. modelsummary #
modelsummary::modelsummary(
    list("EconR Logit" = econ_compat),
    output = "data.frame"
)



#### Dev Tool Testing ####
devtools::load_all()
devtools::test()

devtools::document()

devtools::check()
