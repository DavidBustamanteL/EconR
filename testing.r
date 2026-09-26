################################# Testing Own Econ Package #################################
############################################################################################

# let us begin as usual by setting our wd
rm(list = ls())
gc()

setwd("G:/My Drive/5. CURSOS ONLINE/R Practice FULL/Econ Own Package")

# English Format and Avoiding scientific notation
Sys.setlocale("LC_TIME", "English")
options(scipen = 999)

# Packages
library(tidyverse)
library(magrittr)
library(estimatr)
library(sandwich)
library(clubSandwich)

# Installing test pack
install.packages(
    "G:/My Drive/5. CURSOS ONLINE/R Practice FULL/Econ Own Package",
    repos = NULL,
    type = "source"
)

library(EconR)


#### 1. Basic OLS ####

## With mtcars ##

# EconR
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars)

summary(m_econ)

# estimatr
m_econ_comp = lm_robust(
    mpg ~ wt + hp,
    data = mtcars)

summary(m_econ_comp)

# Coefs comparisson
all.equal(
coef(m_econ),
coef(m_econ_comp))


## testing HC ##

# EconR
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    vcov = "HC3"
)

summary(m_econ)

# estimatr
m_econ_comp = lm_robust(
    mpg ~ wt + hp,
    data = mtcars,
    se_type = "HC3")

summary(m_econ_comp)

# Compare standard errors
sqrt(diag(vcov(m_econ)))
sqrt(diag(vcov(m_econ_comp)))

# Compare covariance matrices
all.equal(
    vcov(m_econ),
    vcov(m_econ_comp)
)

# Compare coefficients
all.equal(
    coef(m_econ),
    coef(m_econ_comp)
)

cbind(
    EconR_coef = coef(m_econ),
    estimatr_coef = coef(m_econ_comp),
    EconR_SE = sqrt(diag(vcov(m_econ))),
    estimatr_SE = sqrt(diag(vcov(m_econ_comp)))
)


## NAs Test ##
df = data.frame(
    y  = c(10, 13, 14, 19, 18, 25, 21, 29, 27, 32),
    x1 = c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10),
    x2 = c(8, 15, 31, 27, NA, 45, 39, 62, 51, 73),
    x3 = c(2, 5, 4, NA, 6, 8, 3, 7, 9, 11)
)

# EconR
m1 = econ_lm(
    y ~ x1 + x2 + x3,
    data = df,
    vcov = "HC3"
)

# estimatr
m2 = lm_robust(
    y ~ x1 + x2 + x3,
    data = df,
    se_type = "HC3"
)

# comparisson
all.equal(
    coef(m1),
    coef(m2)
)

all.equal(
    vcov(m1),
    vcov(m2)
)

cbind(
    EconR_coef = coef(m1),
    estimatr_coef = coef(m2),
    EconR_SE = sqrt(diag(vcov(m1))),
    estimatr_SE = sqrt(diag(vcov(m2)))
)


## Interactions Test ##
m_econ = econ_lm(
    mpg ~ wt * hp,
    data = mtcars,
    vcov = "HC3"
)

m_lm = lm_robust(
    mpg ~ wt * hp,
    data = mtcars,
    se_type = "HC3"
)

# comparisson
all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    vcov(m_econ),
    vcov(m_lm)
)

cbind(
    EconR_coef = coef(m_econ),
    estimatr_coef = coef(m_lm),
    EconR_SE = sqrt(diag(vcov(m_econ))),
    estimatr_SE = sqrt(diag(vcov(m_lm)))
)


## Factor test ##

# Create factor
mtcars$cyl_f = factor(mtcars$cyl)

# EconR
m_econ = econ_lm(
    mpg ~ wt + cyl_f,
    data = mtcars,
    vcov = "HC3"
)

# estimatr
m_comp = lm_robust(
    mpg ~ wt + cyl_f,
    data = mtcars,
    se_type = "HC3"
)

# Comparisson
all.equal(
    coef(m_econ),
    coef(m_comp)
)

all.equal(
    vcov(m_econ),
    vcov(m_comp)
)

cbind(
    EconR_coef = coef(m_econ),
    estimatr_coef = coef(m_comp),
    EconR_SE = sqrt(diag(vcov(m_econ))),
    estimatr_SE = sqrt(diag(vcov(m_comp)))
)


#### 2. WLS ####

## Testing weights ##

# Create weights
mtcars %<>%
    mutate(weighting = mtcars$disp / mean(mtcars$disp))

# EconR
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    weights = weighting,
    vcov = "HC3"
)

# estimatr
m_comp = lm_robust(
    mpg ~ wt + hp,
    data = mtcars,
    weights = weighting,
    se_type = "HC3"
)

# Compare coefficients
all.equal(
    coef(m_econ),
    coef(m_comp)
)

# Compare covariance matrices
all.equal(
    vcov(m_econ),
    vcov(m_comp)
)

# Compare results
cbind(
    EconR_coef = coef(m_econ),
    estimatr_coef = coef(m_comp),
    EconR_SE = sqrt(diag(vcov(m_econ))),
    estimatr_SE = sqrt(diag(vcov(m_comp)))
)


## Testing NAs ##

# adding NAs
df_wls = mtcars %>%
    mutate(
        weighting = disp / mean(disp),
        wt = replace(wt, c(3, 10), NA),
        hp = replace(hp, 5, NA),
        weighting = replace(weighting, 7, NA)
    )

# EconR 
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_wls,
    weights = weighting,
    vcov = "HC3"
)

# estimatr
m_comp = lm_robust(
    mpg ~ wt + hp,
    data = df_wls,
    weights = weighting,
    se_type = "HC3"
)

# Comparing number of observations
c(
    EconR = nobs(m_econ),
    estimatr = nobs(m_comp)
)

# Comparing coefficients
all.equal(
    coef(m_econ),
    coef(m_comp)
)

# Comparing covariance matrices
all.equal(
    vcov(m_econ),
    vcov(m_comp)
)

# Comparing coefficients and standard errors
cbind(
    EconR_coef = coef(m_econ),
    estimatr_coef = coef(m_comp),
    EconR_SE = sqrt(diag(vcov(m_econ))),
    estimatr_SE = sqrt(diag(vcov(m_comp)))
)


#### 3. Final Fit Tests ####

# Fit EconR models
econ_HC0 = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    vcov = "HC0"
)

econ_HC1 = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    vcov = "HC1"
)

econ_HC2 = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    vcov = "HC2"
)

econ_HC3 = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    vcov = "HC3"
)

# Fit estimatr models
comp_HC0 = lm_robust(
    mpg ~ wt + hp,
    data = mtcars,
    se_type = "HC0"
)

comp_HC1 = lm_robust(
    mpg ~ wt + hp,
    data = mtcars,
    se_type = "HC1"
)

comp_HC2 = lm_robust(
    mpg ~ wt + hp,
    data = mtcars,
    se_type = "HC2"
)

comp_HC3 = lm_robust(
    mpg ~ wt + hp,
    data = mtcars,
    se_type = "HC3"
)

# Comparing covariance matrices
c(
    HC0 = isTRUE(all.equal(
        vcov(econ_HC0),
        vcov(comp_HC0)
    )),
    HC1 = isTRUE(all.equal(
        vcov(econ_HC1),
        vcov(comp_HC1)
    )),
    HC2 = isTRUE(all.equal(
        vcov(econ_HC2),
        vcov(comp_HC2)
    )),
    HC3 = isTRUE(all.equal(
        vcov(econ_HC3),
        vcov(comp_HC3)
    ))
)

# Comparing coefficients
c(
    HC0 = isTRUE(all.equal(
        coef(econ_HC0),
        coef(comp_HC0)
    )),
    HC1 = isTRUE(all.equal(
        coef(econ_HC1),
        coef(comp_HC1)
    )),
    HC2 = isTRUE(all.equal(
        coef(econ_HC2),
        coef(comp_HC2)
    )),
    HC3 = isTRUE(all.equal(
        coef(econ_HC3),
        coef(comp_HC3)
    ))
)


#### 4. Clustering Tests ####

# creating cluster var
df_cluster = mtcars %>%
    mutate(
        cluster_id = rep(1:8, each = 4)
    )

## Stata ##

# EconR
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_cluster,
    cluster = cluster_id,
    vcov = "stata"
)

# estimatr
m_comp = lm_robust(
    mpg ~ wt + hp,
    data = df_cluster,
    clusters = cluster_id,
    se_type = "stata"
)

# comp coeffs
all.equal(
    coef(m_econ),
    coef(m_comp)
)

# comp matrix
all.equal(
    vcov(m_econ),
    vcov(m_comp)
)

# comp results
cbind(
    EconR_coef = coef(m_econ),
    estimatr_coef = coef(m_comp),
    EconR_SE = sqrt(diag(vcov(m_econ))),
    estimatr_SE = sqrt(diag(vcov(m_comp)))
)


## CR0 ##

# EconR
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_cluster,
    cluster = cluster_id,
    vcov = "CR0"
)

# estimatr
m_comp = lm_robust(
    mpg ~ wt + hp,
    data = df_cluster,
    clusters = cluster_id,
    se_type = "CR0"
)

# comp coeffs
all.equal(
    coef(m_econ),
    coef(m_comp)
)

# comp matrix
all.equal(
    unname(as.matrix(vcov(m_econ))),
    unname(as.matrix(vcov(m_comp))),
    check.attributes = FALSE
)

# comp results
cbind(
    EconR_coef = coef(m_econ),
    estimatr_coef = coef(m_comp),
    EconR_SE = sqrt(diag(vcov(m_econ))),
    estimatr_SE = sqrt(diag(vcov(m_comp)))
)

## CR2 ##

# EconR
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_cluster,
    cluster = cluster_id,
    vcov = "CR2"
)

# estimatr
m_comp = lm_robust(
    mpg ~ wt + hp,
    data = df_cluster,
    clusters = cluster_id,
    se_type = "CR2"
)

# comp coeffs
all.equal(
    coef(m_econ),
    coef(m_comp)
)

# comp matrix
all.equal(
    unname(as.matrix(vcov(m_econ))),
    unname(as.matrix(vcov(m_comp))),
    check.attributes = FALSE
)

# comp results
cbind(
    EconR_coef = coef(m_econ),
    estimatr_coef = coef(m_comp),
    EconR_SE = sqrt(diag(vcov(m_econ))),
    estimatr_SE = sqrt(diag(vcov(m_comp)))
)


#### 5. Clustering with Missing Values and Weights ####

# creating test data
df_cluster_na = mtcars %>%
    mutate(
        cluster_id = rep(1:8, each = 4),
        test_weight = runif(n()),
        wt = replace(wt, c(3, 10), NA),
        hp = replace(hp, 5, NA),
        test_weight = replace(test_weight, 15, NA)
    )


## 5.1 Stata ##

# EconR
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_cluster_na,
    weights = test_weight,
    cluster = cluster_id,
    vcov = "stata"
)

# estimatr
m_comp = lm_robust(
    mpg ~ wt + hp,
    data = df_cluster_na,
    weights = test_weight,
    clusters = cluster_id,
    se_type = "stata"
)

# comp observations
c(
    EconR = nobs(m_econ),
    estimatr = nobs(m_comp)
)

# comp coeffs
all.equal(
    coef(m_econ),
    coef(m_comp)
)

# comp matrix
all.equal(
    unname(as.matrix(vcov(m_econ))),
    unname(as.matrix(vcov(m_comp))),
    check.attributes = FALSE
)


## 5.2 CR0 ##

# EconR
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_cluster_na,
    weights = test_weight,
    cluster = cluster_id,
    vcov = "CR0"
)

# estimatr
m_comp = lm_robust(
    mpg ~ wt + hp,
    data = df_cluster_na,
    weights = test_weight,
    clusters = cluster_id,
    se_type = "CR0"
)

# comp observations
c(
    EconR = nobs(m_econ),
    estimatr = nobs(m_comp)
)

# comp coeffs
all.equal(
    coef(m_econ),
    coef(m_comp)
)

# comp matrix
all.equal(
    unname(as.matrix(vcov(m_econ))),
    unname(as.matrix(vcov(m_comp))),
    check.attributes = FALSE
)


## 5.3 CR2 ##

# EconR
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_cluster_na,
    weights = test_weight,
    cluster = cluster_id,
    vcov = "CR2"
)

# estimatr
m_comp = lm_robust(
    mpg ~ wt + hp,
    data = df_cluster_na,
    weights = test_weight,
    clusters = cluster_id,
    se_type = "CR2"
)

# comp observations
c(
    EconR = nobs(m_econ),
    estimatr = nobs(m_comp)
)

# comp coeffs
all.equal(
    coef(m_econ),
    coef(m_comp)
)

# comp matrix
all.equal(
    unname(as.matrix(vcov(m_econ))),
    unname(as.matrix(vcov(m_comp))),
    check.attributes = FALSE
)


#### 6. Rank Deficiency Tests ####

# creating perfectly collinear variable
df_rank = mtcars %>%
    mutate(
        wt_copy = wt
    )

## 6.1 Base lm ##
m_lm = lm(
    mpg ~ wt + wt_copy + hp,
    data = df_rank
)
coef(m_lm)
vcov(m_lm)
summary(m_lm)

## 6.2 estimatr ##
m_comp = lm_robust(
    mpg ~ wt + wt_copy + hp,
    data = df_rank,
    se_type = "HC3"
)
coef(m_comp)
vcov(m_comp)
summary(m_comp)

## 6.3 EconR current behavior ##
m_econ = econ_lm(
    mpg ~ wt + wt_copy + hp,
    data = df_rank,
    vcov = "HC3"
)
coef(m_econ)
vcov(m_econ)
summary(m_econ)

## 6.4 Rank Deficiency Across Covariance Types ##
for (type in c("classical", "HC0", "HC1", "HC2", "HC3")) {

    m_econ = econ_lm(
        mpg ~ wt + wt_copy + hp,
        data = df_rank,
        vcov = type
    )

    cat(
        "\n",
        type,
        "\n"
    )

    print(
        coef(m_econ)
    )

    print(
        vcov(m_econ)
    )
}

## 6.5 Rank Deficiency with Clustering ##
df_rank_cluster = df_rank %>%
    mutate(
        cluster_id = rep(1:8, each = 4)
    )

for (type in c("CR0", "CR2", "stata")) {

    # EconR
    m_econ = econ_lm(
        mpg ~ wt + wt_copy + hp,
        data = df_rank_cluster,
        cluster = cluster_id,
        vcov = type
    )

    # estimatr
    m_comp = lm_robust(
        mpg ~ wt + wt_copy + hp,
        data = df_rank_cluster,
        clusters = cluster_id,
        se_type = type
    )

    cat(
        "\n",
        type,
        "\n"
    )

    # coefficients
    print(
        all.equal(
            coef(m_econ),
            coef(m_comp)
        )
    )

    # covariance matrix
    print(
        all.equal(
            unname(as.matrix(vcov(m_econ))),
            unname(as.matrix(vcov(m_comp))),
            check.attributes = FALSE
        )
    )

    # standard errors
    print(
        cbind(
            EconR = sqrt(diag(vcov(m_econ))),
            estimatr = sqrt(diag(vcov(m_comp)))
        )
    )
}


#### 7. No-Intercept Models ####

# Basic no-intercept model
m_econ = econ_lm(
    mpg ~ wt + hp - 1,
    data = mtcars,
    vcov = "HC1"
)

m_ref = estimatr::lm_robust(
    mpg ~ wt + hp - 1,
    data = mtcars,
    se_type = "HC1"
)

coef(m_econ)
coef(m_ref)

all.equal(
    unname(coef(m_econ)),
    unname(coef(m_ref))
)

all.equal(
    unname(vcov(m_econ)),
    unname(vcov(m_ref))
)

# Alternative no-intercept syntax

m_zero = econ_lm(
    mpg ~ 0 + wt + hp,
    data = mtcars,
    vcov = "HC1"
)

all.equal(
    coef(m_econ),
    coef(m_zero)
)

broom::tidy(m_econ)
broom::glance(m_econ)
broom::glance(m_ref)


#### 8. Factors and Contrasts ####

# Factor predictor
df_factor = mtcars %>%
    mutate(
        cyl_factor = factor(cyl)
    )

m_econ = econ_lm(
    mpg ~ cyl_factor + wt,
    data = df_factor,
    vcov = "HC1"
)

m_ref = estimatr::lm_robust(
    mpg ~ cyl_factor + wt,
    data = df_factor,
    se_type = "HC1"
)

coef(m_econ)
coef(m_ref)

all.equal(
    unname(coef(m_econ)),
    unname(coef(m_ref))
)

all.equal(
    unname(vcov(m_econ)),
    unname(vcov(m_ref))
)

model.matrix(m_econ)
model.matrix(
    mpg ~ cyl_factor + wt,
    data = df_factor
)

broom::tidy(m_econ)
broom::tidy(m_ref)


# Different reference level
df_factor_ref = df_factor %>%
    mutate(
        cyl_factor = relevel(
            cyl_factor,
            ref = "8"
        )
    )

m_econ_ref = econ_lm(
    mpg ~ cyl_factor + wt,
    data = df_factor_ref,
    vcov = "HC1"
)

m_ref_ref = estimatr::lm_robust(
    mpg ~ cyl_factor + wt,
    data = df_factor_ref,
    se_type = "HC1"
)

coef(m_econ_ref)
coef(m_ref_ref)

all.equal(
    unname(coef(m_econ_ref)),
    unname(coef(m_ref_ref))
)

all.equal(
    unname(vcov(m_econ_ref)),
    unname(vcov(m_ref_ref))
)


# Prediction with factor levels
new_factor = data.frame(
    cyl_factor = factor(
        c(4, 6, 8),
        levels = levels(df_factor$cyl_factor)
    ),
    wt = c(2.5, 3.0, 3.5)
)

pred_econ = predict(
    m_econ,
    newdata = new_factor
)

pred_ref = predict(
    lm(
        mpg ~ cyl_factor + wt,
        data = df_factor
    ),
    newdata = new_factor
)

pred_econ
pred_ref

all.equal(
    unname(pred_econ),
    unname(pred_ref)
)

# Unseen factor level -> has to give an error if packages behaves correctly!
new_unseen = data.frame(
    cyl_factor = factor(
        "10"
    ),
    wt = 3
)

predict(
    m_econ,
    newdata = new_unseen
)


#### 9. Formula Transformations ####

# Logarithmic and polynomial transformations
m_econ = econ_lm(
    mpg ~ log(wt) + I(hp^2),
    data = mtcars,
    vcov = "HC1"
)

m_ref = estimatr::lm_robust(
    mpg ~ log(wt) + I(hp^2),
    data = mtcars,
    se_type = "HC1"
)

coef(m_econ)
coef(m_ref)

all.equal(
    unname(coef(m_econ)),
    unname(coef(m_ref))
)

all.equal(
    unname(vcov(m_econ)),
    unname(vcov(m_ref))
)

all.equal(
    unname(fitted(m_econ)),
    unname(fitted(m_ref))
)

broom::tidy(m_econ)
broom::tidy(m_ref)


# Prediction with formula transformations
m_lm = lm(
    mpg ~ log(wt) + I(hp^2),
    data = mtcars
)

new_transform = data.frame(
    wt = c(2.5, 3.0, 3.5),
    hp = c(100, 150, 200)
)

pred_econ = predict(
    m_econ,
    newdata = new_transform
)

pred_ref = predict(
    m_lm,
    newdata = new_transform
)

pred_econ
pred_ref

all.equal(
    unname(pred_econ),
    unname(pred_ref)
)


# Invalid value in formula transformation -> Warning NANs if code correct!
df_transform = mtcars %>%
    mutate(
        wt = replace(wt, 1, -1)
    )

m_transform = econ_lm(
    mpg ~ log(wt) + hp,
    data = df_transform
)

nobs(m_transform)

m_transform_ref = lm(
    mpg ~ log(wt) + hp,
    data = df_transform
)

nobs(m_transform)
nobs(m_transform_ref)

coef(m_transform)
coef(m_transform_ref)

all.equal(
    unname(coef(m_transform)),
    unname(coef(m_transform_ref))
)


#### 10. Subset, Missingness, and Weights ####

set.seed(123)

df_combo = mtcars %>%
    mutate(
        test_weight = runif(n()),
        unused = NA_real_,
        wt = replace(wt, c(3, 10), NA),
        hp = replace(hp, 5, NA),
        test_weight = replace(test_weight, 15, NA)
    )

# EconR model
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_combo,
    weights = test_weight,
    subset = cyl != 4,
    vcov = "HC1"
)

# Reference model
m_ref = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_combo,
    weights = test_weight,
    subset = cyl != 4,
    se_type = "HC1"
)

# Number of observations
nobs(m_econ)
nobs(m_ref)

# Coefficients
coef(m_econ)
coef(m_ref)

all.equal(
    unname(coef(m_econ)),
    unname(coef(m_ref))
)

# HC1 covariance matrix
all.equal(
    unname(vcov(m_econ)),
    unname(vcov(m_ref))
)

# Fitted values
all.equal(
    unname(fitted(m_econ)),
    unname(fitted(m_ref))
)


#### 11. Weight Edge Cases ####

# Zero weight
df_weight = mtcars %>%
    mutate(
        test_weight = 1,
        test_weight = replace(
            test_weight,
            1,
            0
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_weight,
    weights = test_weight,
    vcov = "classical"
)

m_ref = lm(
    mpg ~ wt + hp,
    data = df_weight,
    weights = test_weight
)

nobs(m_econ)
nobs(m_ref)

coef(m_econ)
coef(m_ref)

all.equal(
    unname(coef(m_econ)),
    unname(coef(m_ref))
)

fitted(m_econ)[1]
fitted(m_ref)[1]

residuals(m_econ)[1]
residuals(m_ref)[1]


# Negative weight -> Error expected
df_negative = mtcars %>%
    mutate(
        test_weight = 1,
        test_weight = replace(
            test_weight,
            1,
            -1
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_negative,
    weights = test_weight
)

m_ref = lm(
    mpg ~ wt + hp,
    data = df_negative,
    weights = test_weight
)


# Infinite weight -> Must give an error
df_inf = mtcars %>%
    mutate(
        test_weight = 1,
        test_weight = replace(
            test_weight,
            1,
            Inf
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_inf,
    weights = test_weight
)

m_ref = lm(
    mpg ~ wt + hp,
    data = df_inf,
    weights = test_weight
)

m_robust = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_inf,
    weights = test_weight,
    se_type = "HC1"
)
coef(m_robust)


# NaN weight
df_nan = mtcars %>%
    mutate(
        test_weight = 1,
        test_weight = replace(
            test_weight,
            1,
            NaN
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_nan,
    weights = test_weight,
    vcov = "HC1"
)

m_robust = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_nan,
    weights = test_weight,
    se_type = "HC1"
)

nobs(m_econ)
nobs(m_robust)

coef(m_econ)
coef(m_robust)

all.equal(
    unname(coef(m_econ)),
    unname(coef(m_robust))
)

all.equal(
    unname(vcov(m_econ)),
    unname(vcov(m_robust))
)


#### 12. Cluster Edge Cases ####

# Unequal cluster sizes including singleton clusters
df_cluster_edge = mtcars %>%
    mutate(
        cluster_id = c(
            1,
            2, 2,
            3, 3, 3,
            rep(4, 26)
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_cluster_edge,
    vcov = "CR2",
    cluster = cluster_id
)

m_ref = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_cluster_edge,
    clusters = cluster_id,
    se_type = "CR2"
)

coef(m_econ)
coef(m_ref)

all.equal(
    unname(coef(m_econ)),
    unname(coef(m_ref))
)

all.equal(
    unname(as.matrix(vcov(m_econ))),
    unname(as.matrix(vcov(m_ref))),
    check.attributes = FALSE
)

nobs(m_econ)
nobs(m_ref)


# Missing cluster IDs
df_cluster_na = mtcars %>%
    mutate(
        cluster_id = rep(1:8, each = 4),
        cluster_id = replace(
            cluster_id,
            c(3, 10),
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_cluster_na,
    vcov = "CR2",
    cluster = cluster_id
)

m_ref = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_cluster_na,
    clusters = cluster_id,
    se_type = "CR2"
)

nobs(m_econ)
nobs(m_ref)

coef(m_econ)
coef(m_ref)

all.equal(
    unname(coef(m_econ)),
    unname(coef(m_ref))
)

all.equal(
    unname(as.matrix(vcov(m_econ))),
    unname(as.matrix(vcov(m_ref))),
    check.attributes = FALSE
)


# Only one cluster -> should stop, estimatr is flexible here but stat. wrong
df_one_cluster = mtcars %>%
    mutate(
        cluster_id = 1
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_one_cluster,
    vcov = "CR2",
    cluster = cluster_id
)

m_ref = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_one_cluster,
    clusters = cluster_id,
    se_type = "CR2"
)
coef(m_ref)
vcov(m_ref)
broom::tidy(m_ref)
nobs(m_ref)


# Missing cluster ID on an already-excluded observation
df_cluster_model_na = mtcars %>%
    mutate(
        cluster_id = rep(1:8, each = 4),
        wt = replace(
            wt,
            3,
            NA
        ),
        cluster_id = replace(
            cluster_id,
            3,
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_cluster_model_na,
    vcov = "CR2",
    cluster = cluster_id
)

m_ref = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_cluster_model_na,
    clusters = cluster_id,
    se_type = "CR2"
)

nobs(m_econ)
nobs(m_ref)

coef(m_econ)
coef(m_ref)

all.equal(
    unname(coef(m_econ)),
    unname(coef(m_ref))
)

all.equal(
    unname(as.matrix(vcov(m_econ))),
    unname(as.matrix(vcov(m_ref))),
    check.attributes = FALSE
)


#### 13. Prediction Edge Cases ####

# Prediction from rank-deficient model
df_rank_pred = mtcars %>%
    mutate(
        wt_copy = wt
    )

m_econ = econ_lm(
    mpg ~ wt + wt_copy + hp,
    data = df_rank_pred,
    vcov = "HC1"
)

m_ref = lm(
    mpg ~ wt + wt_copy + hp,
    data = df_rank_pred
)

new_rank = data.frame(
    wt = c(2.5, 3.0, 3.5),
    wt_copy = c(2.5, 3.0, 3.5),
    hp = c(100, 150, 200)
)

predict(
    m_econ,
    newdata = new_rank
)

predict(
    m_ref,
    newdata = new_rank
)

all.equal(
    unname(predict(
        m_econ,
        newdata = new_rank
    )),
    unname(predict(
        m_ref,
        newdata = new_rank
    ))
)


# Missing values in newdata
new_missing = data.frame(
    wt = c(2.5, NA, 3.5),
    hp = c(100, 150, NA)
)

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    vcov = "HC1"
)

m_ref = lm(
    mpg ~ wt + hp,
    data = mtcars
)

predict(
    m_econ,
    newdata = new_missing
)

predict(
    m_ref,
    newdata = new_missing
)


# Extra irrelevant columns in newdata
new_extra = data.frame(
    wt = c(2.5, 3.0, 3.5),
    hp = c(100, 150, 200),
    irrelevant = c("A", "B", "C"),
    all_missing = NA_real_
)

pred_econ = predict(
    m_econ,
    newdata = new_extra
)

pred_ref = predict(
    m_ref,
    newdata = new_extra
)

pred_econ
pred_ref

all.equal(
    unname(pred_econ),
    unname(pred_ref)
)


# Required variable absent from newdata -> errors wanted
new_missing_variable = data.frame(
    wt = c(2.5, 3.0, 3.5)
)

predict(
    m_econ,
    newdata = new_missing_variable
)

predict(
    m_ref,
    newdata = new_missing_variable
)


#### 14. Degenerate Models ####

# Constant outcome
df_constant = mtcars %>%
    mutate(
        y_constant = 5
    )

m_econ = econ_lm(
    y_constant ~ wt + hp,
    data = df_constant,
    vcov = "HC1"
)

m_ref = estimatr::lm_robust(
    y_constant ~ wt + hp,
    data = df_constant,
    se_type = "HC1"
)

coef(m_econ)
coef(m_ref)

vcov(m_econ)
vcov(m_ref)

nobs(m_econ)
nobs(m_ref)

broom::tidy(m_econ)
broom::tidy(m_ref)

generics::glance(m_econ)
broom::glance(m_ref)


# Required predictor entirely missing -> error wanted
df_all_na = mtcars %>%
    mutate(
        bad_predictor = NA_real_
    )

m_econ = econ_lm(
    mpg ~ wt + bad_predictor,
    data = df_all_na,
    vcov = "HC1"
)

m_ref = estimatr::lm_robust(
    mpg ~ wt + bad_predictor,
    data = df_all_na,
    se_type = "HC1"
)
coef(m_ref)
vcov(m_ref)
nobs(m_ref)
broom::tidy(m_ref)


# Only one usable observation
df_one_obs = mtcars %>%
    mutate(
        x_test = NA_real_,
        x_test = replace(
            x_test,
            1,
            1
        )
    )

m_econ = econ_lm(
    mpg ~ x_test,
    data = df_one_obs,
    vcov = "HC1"
)

m_ref = estimatr::lm_robust(
    mpg ~ x_test,
    data = df_one_obs,
    se_type = "HC1"
)

coef(m_econ)
coef(m_ref)

vcov(m_econ)
vcov(m_ref)

m_econ$rank
m_ref$rank

m_econ$df.residual
m_ref$df.residual

nobs(m_econ)
nobs(m_ref)

broom::tidy(m_econ)
broom::tidy(m_ref)

generics::glance(m_econ)
broom::glance(m_ref)


# Exactly enough observations for full rank
df_exact = mtcars %>%
    slice(1:3)

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_exact,
    vcov = "HC1"
)

m_ref = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_exact,
    se_type = "HC1"
)

coef(m_econ)
coef(m_ref)

vcov(m_econ)
vcov(m_ref)

m_econ$rank
m_ref$rank

m_econ$df.residual
m_ref$df.residual

nobs(m_econ)
nobs(m_ref)

broom::tidy(m_econ)
broom::tidy(m_ref)

glance.econ_lm(m_econ)
broom::glance(m_ref)


#### 15. All-Zero Weights ####

df_zero_weights = mtcars %>%
    mutate(
        test_weight = 0
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_zero_weights,
    weights = test_weight,
    vcov = "HC1"
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_zero_weights,
    weights = test_weight
)

m_ref = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_zero_weights,
    weights = test_weight,
    se_type = "HC1"
)

coef(m_lm)
coef(m_ref)

nobs(m_lm)
nobs(m_ref)

m_lm$rank
m_ref$rank

m_lm$df.residual
m_ref$df.residual


#### 16. WLS Goodness-of-Fit ####

df_wls = mtcars %>%
    mutate(
        test_weight = seq(
            0.5,
            2,
            length.out = n()
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_wls,
    weights = test_weight,
    vcov = "HC1"
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_wls,
    weights = test_weight
)

m_ref = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_wls,
    weights = test_weight,
    se_type = "HC1"
)

generics::glance(m_econ)
broom::glance(m_lm)
broom::glance(m_ref)


# WLS Goodness-of-Fit with Zero Weights
df_wls_zero = df_wls %>%
    mutate(
        test_weight = replace(
            test_weight,
            c(1, 5, 10),
            0
        )
    )

m_econ_zero = econ_lm(
    mpg ~ wt + hp,
    data = df_wls_zero,
    weights = test_weight,
    vcov = "HC1"
)

m_lm_zero = lm(
    mpg ~ wt + hp,
    data = df_wls_zero,
    weights = test_weight
)

generics::glance(m_econ_zero)

broom::glance(m_lm_zero)

nobs(m_econ_zero)
nobs(m_lm_zero)

m_econ_zero$df.residual
df.residual(m_lm_zero)

# test inference hc1
m_ref_zero = estimatr::lm_robust(
    mpg ~ wt + hp,
    data = df_wls_zero,
    weights = test_weight,
    se_type = "HC1"
)

coef(m_econ_zero)
coef(m_ref_zero)

sqrt(diag(vcov(m_econ_zero)))
sqrt(diag(vcov(m_ref_zero)))

m_econ_zero$df.residual
m_ref_zero$df.residual


#### 17. NA Exclude Semantics ####

df_na_exclude = mtcars %>%
    mutate(
        wt = replace(
            wt,
            c(3, 8, 15),
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_na_exclude,
    na.action = na.exclude,
    vcov = "HC1"
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_na_exclude,
    na.action = na.exclude
)

nobs(m_econ)
nobs(m_lm)

length(fitted(m_econ))
length(fitted(m_lm))

length(residuals(m_econ))
length(residuals(m_lm))

fitted(m_econ)
fitted(m_lm)

residuals(m_econ)
residuals(m_lm)


#### 18. Offset Terms ####

m_econ = econ_lm(
    mpg ~ wt + hp + offset(disp / 100),
    data = mtcars,
    vcov = "HC1"
)

m_lm = lm(
    mpg ~ wt + hp + offset(disp / 100),
    data = mtcars
)

coef(m_econ)
coef(m_lm)

fitted(m_econ)
fitted(m_lm)

residuals(m_econ)
residuals(m_lm)

# hc1 inference
m_ref = estimatr::lm_robust(
    mpg ~ wt + hp + offset(disp / 100),
    data = mtcars,
    se_type = "HC1"
)

sqrt(diag(vcov(m_econ)))
sqrt(diag(vcov(m_ref)))

all.equal(
    unname(vcov(m_econ)),
    unname(vcov(m_ref)),
    tolerance = 1e-12
)


# new data
new_offset = mtcars %>%
    slice(1:5) %>%
    mutate(
        wt = wt * 1.05,
        disp = disp * 1.10
    )

predict(
    m_econ,
    newdata = new_offset
)

predict(
    m_lm,
    newdata = new_offset
)

all.equal(
    predict(
        m_econ,
        newdata = new_offset
    ),
    predict(
        m_lm,
        newdata = new_offset
    ),
    tolerance = 1e-12
)


#### 19. Custom contracts ####

# non default reference
df_contrast = mtcars %>%
    mutate(
        cyl_factor = factor(cyl),
        cyl_factor = relevel(
            cyl_factor,
            ref = "6"
        )
    )

m_econ = econ_lm(
    mpg ~ wt + cyl_factor,
    data = df_contrast,
    vcov = "HC1"
)

m_lm = lm(
    mpg ~ wt + cyl_factor,
    data = df_contrast
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    model.matrix(m_econ),
    model.matrix(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-12
)

all.equal(
    unname(vcov(m_econ)),
    unname(
        sandwich::vcovHC(
            m_lm,
            type = "HC1"
        )
    ),
    tolerance = 1e-12
)


# explicit sum contrasts
df_contrast_sum = mtcars %>%
    mutate(
        cyl_factor = factor(cyl)
    )

contrasts(df_contrast_sum$cyl_factor) = contr.sum(3)

m_econ = econ_lm(
    mpg ~ wt + cyl_factor,
    data = df_contrast_sum,
    vcov = "HC1"
)

m_lm = lm(
    mpg ~ wt + cyl_factor,
    data = df_contrast_sum
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    model.matrix(m_econ),
    model.matrix(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-12
)

all.equal(
    unname(vcov(m_econ)),
    unname(
        sandwich::vcovHC(
            m_lm,
            type = "HC1"
        )
    ),
    tolerance = 1e-12
)


# prediction with custom contrasts
new_contrast = data.frame(
    wt = c(
        2.5,
        3.0,
        3.5
    ),
    cyl_factor = factor(
        c(
            4,
            6,
            8
        ),
        levels = c(
            4,
            6,
            8
        )
    )
)

predict(
    m_econ,
    newdata = new_contrast
)

predict(
    m_lm,
    newdata = new_contrast
)

all.equal(
    predict(
        m_econ,
        newdata = new_contrast
    ),
    predict(
        m_lm,
        newdata = new_contrast
    ),
    tolerance = 1e-12
)


# ordered factor with polynomial contrasts
df_ordered = mtcars %>%
    mutate(
        cyl_ordered = ordered(
            cyl,
            levels = c(
                4,
                6,
                8
            )
        )
    )

m_econ = econ_lm(
    mpg ~ wt + cyl_ordered,
    data = df_ordered,
    vcov = "HC1"
)

m_lm = lm(
    mpg ~ wt + cyl_ordered,
    data = df_ordered
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    model.matrix(m_econ),
    model.matrix(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-12
)

all.equal(
    unname(vcov(m_econ)),
    unname(
        sandwich::vcovHC(
            m_lm,
            type = "HC1"
        )
    ),
    tolerance = 1e-12
)


# prediction with ordered factor
new_ordered = data.frame(
    wt = c(
        2.5,
        3.0,
        3.5
    ),
    cyl_ordered = ordered(
        c(
            4,
            6,
            8
        ),
        levels = c(
            4,
            6,
            8
        )
    )
)

predict(
    m_econ,
    newdata = new_ordered
)

predict(
    m_lm,
    newdata = new_ordered
)

all.equal(
    predict(
        m_econ,
        newdata = new_ordered
    ),
    predict(
        m_lm,
        newdata = new_ordered
    ),
    tolerance = 1e-12
)


#### 20. Weighted Offset Models ####

df_offset_wls = mtcars %>%
    mutate(
        test_weight = seq(
            0.5,
            2,
            length.out = n()
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp + offset(disp / 100),
    data = df_offset_wls,
    weights = test_weight,
    vcov = "HC1"
)

m_lm = lm(
    mpg ~ wt + hp + offset(disp / 100),
    data = df_offset_wls,
    weights = test_weight
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-12
)

all.equal(
    unname(vcov(m_econ)),
    unname(
        sandwich::vcovHC(
            m_lm,
            type = "HC1"
        )
    ),
    tolerance = 1e-12
)


# weighted offset with missingness
df_offset_missing = mtcars %>%
    mutate(
        test_weight = seq(
            0.5,
            2,
            length.out = n()
        ),
        offset_var = disp / 100,
        test_weight = replace(
            test_weight,
            c(3, 10),
            NA
        ),
        offset_var = replace(
            offset_var,
            c(5, 12),
            NA
        ),
        hp = replace(
            hp,
            c(7, 18),
            NA
        )
    )

m_econ = suppressWarnings(
    econ_lm(
        mpg ~ wt + hp + offset(offset_var),
        data = df_offset_missing,
        weights = test_weight,
        vcov = "HC1"
    )
)

m_lm = lm(
    mpg ~ wt + hp + offset(offset_var),
    data = df_offset_missing,
    weights = test_weight
)

nobs(m_econ)
nobs(m_lm)

coef(m_econ)
coef(m_lm)

all.equal(
    nobs(m_econ),
    nobs(m_lm)
)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-12
)

all.equal(
    unname(vcov(m_econ)),
    unname(
        sandwich::vcovHC(
            m_lm,
            type = "HC1"
        )
    ),
    tolerance = 1e-12
)


# weighted offset with zero weights
df_offset_zero = mtcars %>%
    mutate(
        test_weight = seq(
            0.5,
            2,
            length.out = n()
        ),
        test_weight = replace(
            test_weight,
            c(1, 5, 10),
            0
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp + offset(disp / 100),
    data = df_offset_zero,
    weights = test_weight,
    vcov = "HC1"
)

m_lm = lm(
    mpg ~ wt + hp + offset(disp / 100),
    data = df_offset_zero,
    weights = test_weight
)

nobs(m_econ)
nobs(m_lm)

m_econ$df.residual
m_lm$df.residual

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-12
)

all.equal(
    unname(vcov(m_econ)),
    unname(
        sandwich::vcovHC(
            m_lm,
            type = "HC1"
        )
    ),
    tolerance = 1e-12
)


#### 21. NA Fail ####

# missing covariate
df_na_fail = mtcars %>%
    mutate(
        hp = replace(
            hp,
            c(3, 8),
            NA
        )
    )

m_econ = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_na_fail,
        na.action = na.fail
    ),
    error = function(e) e
)

m_lm = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_na_fail,
        na.action = na.fail
    ),
    error = function(e) e
)

m_econ
m_lm

conditionMessage(m_econ)
conditionMessage(m_lm)


# missing weights
df_na_fail_weight = mtcars %>%
    mutate(
        test_weight = seq(
            0.5,
            2,
            length.out = n()
        ),
        test_weight = replace(
            test_weight,
            c(3, 8),
            NA
        )
    )

m_econ = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_na_fail_weight,
        weights = test_weight,
        na.action = na.fail
    ),
    error = function(e) e
)

m_lm = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_na_fail_weight,
        weights = test_weight,
        na.action = na.fail
    ),
    error = function(e) e
)

m_econ
m_lm

conditionMessage(m_econ)
conditionMessage(m_lm)


# missing cluster IDs
df_na_fail_cluster = mtcars %>%
    mutate(
        cluster_id = rep(
            1:8,
            each = 4
        ),
        cluster_id = replace(
            cluster_id,
            c(3, 8),
            NA
        )
    )

m_econ = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_na_fail_cluster,
        cluster = cluster_id,
        vcov = "CR0",
        na.action = na.fail
    ),
    error = function(e) e
)

m_econ
conditionMessage(m_econ)


# na fail
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    na.action = na.fail
)

m_lm = lm(
    mpg ~ wt + hp,
    data = mtcars,
    na.action = na.fail
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)


#### 22. Subset edge cases ####

# Basic Subset
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = cyl == 6
)

m_lm = lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = cyl == 6
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)

rownames(model.frame(m_econ))
rownames(model.frame(m_lm))


# Subset with missing values outside estimation sample
df_subset_na = mtcars %>%
    mutate(
        hp = replace(
            hp,
            which(cyl != 6)[1:3],
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_subset_na,
    subset = cyl == 6
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_subset_na,
    subset = cyl == 6
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)


# Subset with missing value inside estimation sample
df_subset_na_inside = mtcars %>%
    mutate(
        hp = replace(
            hp,
            which(cyl == 6)[1],
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_subset_na_inside,
    subset = cyl == 6
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_subset_na_inside,
    subset = cyl == 6
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)

rownames(model.frame(m_econ))
rownames(model.frame(m_lm))


# Subset with missing weights outside estimation sample
df_subset_weight = mtcars %>%
    mutate(
        test_weight = seq(
            0.5,
            2,
            length.out = n()
        ),
        test_weight = replace(
            test_weight,
            which(cyl != 6)[1],
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_subset_weight,
    subset = cyl == 6,
    weights = test_weight
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_subset_weight,
    subset = cyl == 6,
    weights = test_weight
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)


# Subset with missing weight inside estimation sample -> warning expected
df_subset_weight_inside = mtcars %>%
    mutate(
        test_weight = seq(
            0.5,
            2,
            length.out = n()
        ),
        test_weight = replace(
            test_weight,
            which(cyl == 6)[1],
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_subset_weight_inside,
    subset = cyl == 6,
    weights = test_weight
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_subset_weight_inside,
    subset = cyl == 6,
    weights = test_weight
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)

rownames(model.frame(m_econ))
rownames(model.frame(m_lm))


# Subset with missing cluster ID outside estimation sample
df_subset_cluster = mtcars %>%
    mutate(
        cluster_id = rep(
            1:8,
            each = 4
        ),
        cluster_id = replace(
            cluster_id,
            which(cyl != 6)[1],
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_subset_cluster,
    subset = cyl == 6,
    cluster = cluster_id,
    vcov = "CR0"
)

m_full = econ_lm(
    mpg ~ wt + hp,
    data = df_subset_cluster %>%
        filter(cyl == 6),
    cluster = cluster_id,
    vcov = "CR0"
)

all.equal(
    coef(m_econ),
    coef(m_full)
)

all.equal(
    vcov(m_econ),
    vcov(m_full)
)

nobs(m_econ)
nobs(m_full)


# Subset with missing cluster ID inside estimation sample
df_subset_cluster_inside = mtcars %>%
    mutate(
        cluster_id = rep(
            1:8,
            each = 4
        ),
        cluster_id = replace(
            cluster_id,
            which(cyl == 6)[1],
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_subset_cluster_inside,
    subset = cyl == 6,
    cluster = cluster_id,
    vcov = "CR0"
)

m_full = econ_lm(
    mpg ~ wt + hp,
    data = df_subset_cluster_inside %>%
        filter(
            cyl == 6,
            !is.na(cluster_id)
        ),
    cluster = cluster_id,
    vcov = "CR0"
)

all.equal(
    coef(m_econ),
    coef(m_full)
)

all.equal(
    vcov(m_econ),
    vcov(m_full)
)

nobs(m_econ)
nobs(m_full)

rownames(model.frame(m_econ))
rownames(model.frame(m_full))


# na.fail with missing value outside subset
df_subset_fail = mtcars %>%
    mutate(
        hp = replace(
            hp,
            which(cyl != 6)[1],
            NA
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_subset_fail,
    subset = cyl == 6,
    na.action = na.fail
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_subset_fail,
    subset = cyl == 6,
    na.action = na.fail
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

nobs(m_econ)
nobs(m_lm)


# na.fail with missing value inside subset -> warning expected but in plain text
df_subset_fail_inside = mtcars %>%
    mutate(
        hp = replace(
            hp,
            which(cyl == 6)[1],
            NA
        )
    )

m_econ = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_subset_fail_inside,
        subset = cyl == 6,
        na.action = na.fail
    ),
    error = function(e) e
)

m_lm = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_subset_fail_inside,
        subset = cyl == 6,
        na.action = na.fail
    ),
    error = function(e) e
)

conditionMessage(m_econ)
conditionMessage(m_lm)


#### 23. Formula and Environment Edge Cases ####

# Predictor defined outside data
external_x = mtcars$wt * 2

m_econ = econ_lm(
    mpg ~ external_x + hp,
    data = mtcars
)

m_lm = lm(
    mpg ~ external_x + hp,
    data = mtcars
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)

rownames(model.frame(m_econ))
rownames(model.frame(m_lm))


# Custom function defined outside data
double_value = function(x) {
    x * 2
}

m_econ = econ_lm(
    mpg ~ double_value(wt) + hp,
    data = mtcars
)

m_lm = lm(
    mpg ~ double_value(wt) + hp,
    data = mtcars
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    model.matrix(m_econ),
    model.matrix(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)


# External scalar used in formula
scale_factor = 3.5

m_econ = econ_lm(
    mpg ~ I(wt * scale_factor) + hp,
    data = mtcars
)

m_lm = lm(
    mpg ~ I(wt * scale_factor) + hp,
    data = mtcars
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    model.matrix(m_econ),
    model.matrix(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)


# Weights defined outside data
external_weights = seq(
    0.5,
    2,
    length.out = nrow(mtcars)
)

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    weights = external_weights
)

m_lm = lm(
    mpg ~ wt + hp,
    data = mtcars,
    weights = external_weights
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

all.equal(
    weights(m_econ),
    weights(m_lm)
)

nobs(m_econ)
nobs(m_lm)


# Subset defined outside data
external_subset = mtcars$cyl == 6

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = external_subset
)

m_lm = lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = external_subset
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

nobs(m_econ)
nobs(m_lm)

rownames(model.frame(m_econ))
rownames(model.frame(m_lm))


# Cluster variable defined outside data
external_cluster = rep(
    1:8,
    each = 4
)

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    cluster = external_cluster,
    vcov = "CR0"
)

m_internal = econ_lm(
    mpg ~ wt + hp,
    data = mtcars %>%
        mutate(
            cluster_id = external_cluster
        ),
    cluster = cluster_id,
    vcov = "CR0"
)

all.equal(
    coef(m_econ),
    coef(m_internal)
)

all.equal(
    vcov(m_econ),
    vcov(m_internal)
)

all.equal(
    m_econ$cluster,
    m_internal$cluster
)

nobs(m_econ)
nobs(m_internal)


# Locally scoped formula environment
fit_local_models = function() {

    local_multiplier = 2.5

    local_transform = function(x) {
        x * local_multiplier
    }

    m_econ = econ_lm(
        mpg ~ local_transform(wt) + hp,
        data = mtcars
    )

    m_lm = lm(
        mpg ~ local_transform(wt) + hp,
        data = mtcars
    )

    list(
        econ = m_econ,
        lm = m_lm
    )
}

local_models = fit_local_models()

all.equal(
    coef(local_models$econ),
    coef(local_models$lm)
)

all.equal(
    model.matrix(local_models$econ),
    model.matrix(local_models$lm)
)

all.equal(
    fitted(local_models$econ),
    fitted(local_models$lm)
)

all.equal(
    residuals(local_models$econ),
    residuals(local_models$lm)
)

nobs(local_models$econ)
nobs(local_models$lm)


#### 24. Prediction Edge Cases ####

# Missing predictor values in newdata
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = mtcars
)

m_lm = lm(
    mpg ~ wt + hp,
    data = mtcars
)

newdata_na = data.frame(
    wt = c(
        2.5,
        NA,
        3.5
    ),
    hp = c(
        100,
        120,
        NA
    )
)

pred_econ = predict(
    m_econ,
    newdata = newdata_na
)

pred_lm = predict(
    m_lm,
    newdata = newdata_na
)

pred_econ
pred_lm

all.equal(
    pred_econ,
    pred_lm
)

is.na(pred_econ)
is.na(pred_lm)


# Missing predictor column in newdata -> Plain text warning wanted
newdata_missing = data.frame(
    wt = c(
        2.5,
        3.0,
        3.5
    )
)

error_econ = tryCatch(
    predict(
        m_econ,
        newdata = newdata_missing
    ),
    error = function(e) e
)

error_lm = tryCatch(
    predict(
        m_lm,
        newdata = newdata_missing
    ),
    error = function(e) e
)

conditionMessage(error_econ)
conditionMessage(error_lm)


# Unseen factor level in newdata -> Plain text warning wanted
df_factor = mtcars %>%
    mutate(
        cyl_factor = factor(cyl)
    )

m_econ_factor = econ_lm(
    mpg ~ wt + cyl_factor,
    data = df_factor
)

m_lm_factor = lm(
    mpg ~ wt + cyl_factor,
    data = df_factor
)

newdata_new_level = data.frame(
    wt = 3,
    cyl_factor = factor(
        10,
        levels = c(
            4,
            6,
            8,
            10
        )
    )
)

error_econ = tryCatch(
    predict(
        m_econ_factor,
        newdata = newdata_new_level
    ),
    error = function(e) e
)

error_lm = tryCatch(
    predict(
        m_lm_factor,
        newdata = newdata_new_level
    ),
    error = function(e) e
)

conditionMessage(error_econ)
conditionMessage(error_lm)


# Valid factor levels in newdata
newdata_factor = data.frame(
    wt = c(
        2.5,
        3.0,
        3.5
    ),
    cyl_factor = factor(
        c(
            4,
            6,
            8
        ),
        levels = c(
            4,
            6,
            8
        )
    )
)

pred_econ = predict(
    m_econ_factor,
    newdata = newdata_factor
)

pred_lm = predict(
    m_lm_factor,
    newdata = newdata_factor
)

pred_econ
pred_lm

all.equal(
    pred_econ,
    pred_lm
)


# Prediction with transformed terms
m_econ_transform = econ_lm(
    mpg ~ log(wt) + I(hp^2),
    data = mtcars
)

m_lm_transform = lm(
    mpg ~ log(wt) + I(hp^2),
    data = mtcars
)

newdata_transform = data.frame(
    wt = c(
        2.5,
        3.0,
        3.5
    ),
    hp = c(
        100,
        150,
        200
    )
)

pred_econ = predict(
    m_econ_transform,
    newdata = newdata_transform
)

pred_lm = predict(
    m_lm_transform,
    newdata = newdata_transform
)

pred_econ
pred_lm

all.equal(
    pred_econ,
    pred_lm
)


# Prediction after rank deficiency
df_rank_pred = mtcars %>%
    mutate(
        wt_copy = wt
    )

m_econ_rank = econ_lm(
    mpg ~ wt + wt_copy + hp,
    data = df_rank_pred
)

m_lm_rank = lm(
    mpg ~ wt + wt_copy + hp,
    data = df_rank_pred
)

newdata_rank = data.frame(
    wt = c(
        2.5,
        3.0,
        3.5
    ),
    wt_copy = c(
        2.5,
        3.0,
        3.5
    ),
    hp = c(
        100,
        150,
        200
    )
)

pred_econ = predict(
    m_econ_rank,
    newdata = newdata_rank
)

pred_lm = predict(
    m_lm_rank,
    newdata = newdata_rank
)

pred_econ
pred_lm

all.equal(
    pred_econ,
    pred_lm
)


# Prediction with interaction terms
m_econ_interaction = econ_lm(
    mpg ~ wt * hp,
    data = mtcars
)

m_lm_interaction = lm(
    mpg ~ wt * hp,
    data = mtcars
)

newdata_interaction = data.frame(
    wt = c(
        2.5,
        3.0,
        3.5
    ),
    hp = c(
        100,
        150,
        200
    )
)

pred_econ = predict(
    m_econ_interaction,
    newdata = newdata_interaction
)

pred_lm = predict(
    m_lm_interaction,
    newdata = newdata_interaction
)

pred_econ
pred_lm

all.equal(
    pred_econ,
    pred_lm
)


# Prediction with na.omit
newdata_omit = data.frame(
    wt = c(
        2.5,
        NA,
        3.5
    ),
    hp = c(
        100,
        120,
        200
    )
)

pred_econ = predict(
    m_econ,
    newdata = newdata_omit,
    na.action = na.omit
)

pred_lm = predict(
    m_lm,
    newdata = newdata_omit,
    na.action = na.omit
)

pred_econ
pred_lm

all.equal(
    pred_econ,
    pred_lm
)

length(pred_econ)
length(pred_lm)

names(pred_econ)
names(pred_lm)


#### 25. Degrees of Freedom and Tiny Samples ####

# Exactly saturated model
df_saturated = data.frame(
    y = c(
        2,
        5,
        9
    ),
    x1 = c(
        1,
        2,
        4
    ),
    x2 = c(
        2,
        5,
        3
    )
)

m_econ = econ_lm(
    y ~ x1 + x2,
    data = df_saturated
)

m_lm = lm(
    y ~ x1 + x2,
    data = df_saturated
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

df.residual(m_econ)
df.residual(m_lm)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)


# Summary with zero residual degrees of freedom
summary_econ = summary(m_econ)
summary_lm = summary(m_lm)

summary_econ$coefficients
summary_lm$coefficients

summary_econ$df.residual
summary_lm$df[2]

all.equal(
    summary_econ$coefficients,
    summary_lm$coefficients
)


# One residual degree of freedom
df_df1 = data.frame(
    y = c(
        2,
        5,
        9,
        10
    ),
    x1 = c(
        1,
        2,
        4,
        5
    ),
    x2 = c(
        2,
        5,
        3,
        7
    )
)

m_econ = econ_lm(
    y ~ x1 + x2,
    data = df_df1
)

m_lm = lm(
    y ~ x1 + x2,
    data = df_df1
)

summary_econ = summary(m_econ)
summary_lm = summary(m_lm)

df.residual(m_econ)
df.residual(m_lm)

summary_econ$coefficients
summary_lm$coefficients

all.equal(
    summary_econ$coefficients,
    summary_lm$coefficients
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)


# Single-observation intercept-only model
df_single = data.frame(
    y = 5
)

m_econ = econ_lm(
    y ~ 1,
    data = df_single
)

m_lm = lm(
    y ~ 1,
    data = df_single
)

coef(m_econ)
coef(m_lm)

df.residual(m_econ)
df.residual(m_lm)

fitted(m_econ)
fitted(m_lm)

residuals(m_econ)
residuals(m_lm)

summary(m_econ)$coefficients
summary(m_lm)$coefficients


# Constant outcome
df_constant_y = mtcars %>%
    mutate(
        constant_y = 10
    )

m_econ = econ_lm(
    constant_y ~ wt + hp,
    data = df_constant_y
)

m_lm = lm(
    constant_y ~ wt + hp,
    data = df_constant_y
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

glance(m_econ)

summary(m_lm)$r.squared
summary(m_lm)$adj.r.squared

# Constant-outcome diagnostic
y_econ = model.response(
    model.frame(m_econ)
)

resid_econ = residuals(
    m_econ
)

fitted_econ = fitted(
    m_econ
)

y_lm = model.response(
    model.frame(m_lm)
)

resid_lm = residuals(
    m_lm
)

fitted_lm = fitted(
    m_lm
)

c(
    econ_rss = sum(resid_econ^2),
    lm_rss = sum(resid_lm^2),
    econ_tss = sum((y_econ - mean(y_econ))^2),
    lm_tss = sum((y_lm - mean(y_lm))^2),
    econ_fitted_variation = sum(
        (fitted_econ - mean(fitted_econ))^2
    ),
    lm_fitted_variation = sum(
        (fitted_lm - mean(fitted_lm))^2
    )
)

range(resid_econ)
range(resid_lm)

range(fitted_econ)
range(fitted_lm)


# Constant predictor
df_constant_x = mtcars %>%
    mutate(
        constant_x = 5
    )

m_econ = econ_lm(
    mpg ~ constant_x + wt,
    data = df_constant_x
)

m_lm = lm(
    mpg ~ constant_x + wt,
    data = df_constant_x
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

m_econ$rank
m_lm$rank

df.residual(m_econ)
df.residual(m_lm)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)


# Intercept-only model
m_econ = econ_lm(
    mpg ~ 1,
    data = mtcars
)

m_lm = lm(
    mpg ~ 1,
    data = mtcars
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

df.residual(m_econ)
df.residual(m_lm)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

summary(m_econ)$coefficients
summary(m_lm)$coefficients


# No-intercept model
m_econ = econ_lm(
    mpg ~ 0 + wt + hp,
    data = mtcars
)

m_lm = lm(
    mpg ~ 0 + wt + hp,
    data = mtcars
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

df.residual(m_econ)
df.residual(m_lm)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    residuals(m_econ),
    residuals(m_lm)
)

summary(m_econ)$coefficients
summary(m_lm)$coefficients

glance(m_econ)

summary(m_lm)$r.squared
summary(m_lm)$adj.r.squared


# saturated model
df_saturated = data.frame(
    y = c(2, 5, 9),
    x1 = c(1, 2, 4),
    x2 = c(2, 5, 3)
)

m_econ = econ_lm(
    y ~ x1 + x2,
    data = df_saturated
)

m_lm = lm(
    y ~ x1 + x2,
    data = df_saturated
)

tidy_econ = tidy(
    m_econ
)

tidy_lm = broom::tidy(
    m_lm,
    conf.int = TRUE
)

tidy_econ
tidy_lm


# glance() with saturated model
glance_econ = glance(
    m_econ
)

glance_econ

c(
    r_squared = summary(m_lm)$r.squared,
    adj_r_squared = summary(m_lm)$adj.r.squared,
    aic = AIC(m_lm),
    bic = BIC(m_lm),
    nobs = nobs(m_lm)
)


# WLS with one residual DoF
df_df1_wls = data.frame(
    y = c(2, 5, 9, 10),
    x1 = c(1, 2, 4, 5),
    x2 = c(2, 5, 3, 7),
    w = c(1, 2, 3, 4)
)

m_econ = econ_lm(
    y ~ x1 + x2,
    data = df_df1_wls,
    weights = w
)

m_lm = lm(
    y ~ x1 + x2,
    data = df_df1_wls,
    weights = w
)

df.residual(m_econ)
df.residual(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

summary(m_econ)$coefficients
summary(m_lm)$coefficients

all.equal(
    summary(m_econ)$coefficients,
    summary(m_lm)$coefficients
)


# WLS saturated model (df.residual = 0)
df_saturated_wls = data.frame(
    y = c(2, 5, 9),
    x1 = c(1, 2, 4),
    x2 = c(2, 5, 3),
    w = c(1, 2, 3)
)

m_econ = econ_lm(
    y ~ x1 + x2,
    data = df_saturated_wls,
    weights = w
)

m_lm = lm(
    y ~ x1 + x2,
    data = df_saturated_wls,
    weights = w
)

df.residual(m_econ)
df.residual(m_lm)

coef(m_econ)
coef(m_lm)

summary(m_econ)$coefficients
summary(m_lm)$coefficients

all.equal(
    summary(m_econ)$coefficients,
    summary(m_lm)$coefficients
)


# Degrees of freedom with zero weights
df_zero_w = data.frame(
    y = c(2, 5, 9, 10, 14),
    x = c(1, 2, 4, 5, 7),
    w = c(1, 1, 1, 0, 0)
)

m_econ = econ_lm(
    y ~ x,
    data = df_zero_w,
    weights = w
)

m_lm = lm(
    y ~ x,
    data = df_zero_w,
    weights = w
)

nobs(m_econ)
nobs(m_lm)

df.residual(m_econ)
df.residual(m_lm)

coef(m_econ)
coef(m_lm)

summary(m_econ)$coefficients
summary(m_lm)$coefficients


#### 26. Numerical and Pathological Inputs ####

# Very large predictor values
df_large = mtcars %>%
    dplyr::mutate(
        wt_large = wt * 1e10,
        hp_large = hp * 1e8
    )

m_econ = econ_lm(
    mpg ~ wt_large + hp_large,
    data = df_large
)

m_lm = lm(
    mpg ~ wt_large + hp_large,
    data = df_large
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-12
)

m_econ$rank
m_lm$rank


# Very small predictor values
df_small = mtcars %>%
    dplyr::mutate(
        wt_small = wt * 1e-10,
        hp_small = hp * 1e-8
    )

m_econ = econ_lm(
    mpg ~ wt_small + hp_small,
    data = df_small
)

m_lm = lm(
    mpg ~ wt_small + hp_small,
    data = df_small
)

coef(m_econ)
coef(m_lm)

m_econ$rank
m_lm$rank

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-12
)


# Very large outcome values
df_large_y = mtcars %>%
    dplyr::mutate(
        mpg_large = mpg * 1e12
    )

m_econ = econ_lm(
    mpg_large ~ wt + hp,
    data = df_large_y
)

m_lm = lm(
    mpg_large ~ wt + hp,
    data = df_large_y
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-12
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-12
)


# Very small outcome values
df_small_y = mtcars %>%
    dplyr::mutate(
        mpg_small = mpg * 1e-12
    )

m_econ = econ_lm(
    mpg_small ~ wt + hp,
    data = df_small_y
)

m_lm = lm(
    mpg_small ~ wt + hp,
    data = df_small_y
)

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-12
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-12
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-12
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-12
)


# Extreme but heterogeneous weights
df_extreme_w = mtcars %>%
    dplyr::mutate(
        w_extreme = 10^seq(
            -6,
            6,
            length.out = n()
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_extreme_w,
    weights = w_extreme
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_extreme_w,
    weights = w_extreme
)

coef(m_econ)
coef(m_lm)

m_econ$rank
m_lm$rank

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-10
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-10
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-10
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-10
)


# Near-collinearity
df_near_collinear = mtcars %>%
    dplyr::mutate(
        wt_near = wt + seq_along(wt) * 1e-8
    )

m_econ = econ_lm(
    mpg ~ wt + wt_near + hp,
    data = df_near_collinear
)

m_lm = lm(
    mpg ~ wt + wt_near + hp,
    data = df_near_collinear
)

coef(m_econ)
coef(m_lm)

m_econ$rank
m_lm$rank

m_econ$aliased
which(is.na(coef(m_lm)))

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-10
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-10
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-10
)


# Near-collinearity that remains full rank
df_near_full = mtcars %>%
    dplyr::mutate(
        wt_near = wt + seq_along(wt) * 1e-5
    )

m_econ = econ_lm(
    mpg ~ wt + wt_near + hp,
    data = df_near_full
)

m_lm = lm(
    mpg ~ wt + wt_near + hp,
    data = df_near_full
)

coef(m_econ)
coef(m_lm)

m_econ$rank
m_lm$rank

m_econ$aliased
which(is.na(coef(m_lm)))

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-8
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-8
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-8
)


# Extremely unbalanced weights including zero
df_weight_mix = mtcars %>%
    dplyr::mutate(
        w_mix = c(
            rep(0, 4),
            10^seq(
                -8,
                8,
                length.out = n() - 4
            )
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_weight_mix,
    weights = w_mix
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_weight_mix,
    weights = w_mix
)

nobs(m_econ)
nobs(m_lm)

df.residual(m_econ)
df.residual(m_lm)

m_econ$rank
m_lm$rank

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-8
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-8
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-8
)


# Infinite predictor value -> errors wanted in plain text
df_inf = mtcars %>%
    dplyr::mutate(
        wt = replace(
            wt,
            5,
            Inf
        )
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_inf
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_inf
    ),
    error = function(e) e
)

econ_result
lm_result


# Infinite outcome -> plain text errors wanted
df_inf_y = mtcars %>%
    dplyr::mutate(
        mpg = replace(
            mpg,
            5,
            Inf
        )
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_inf_y
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_inf_y
    ),
    error = function(e) e
)

econ_result
lm_result


# infinite Weight -> same as above
df_inf_w = mtcars %>%
    dplyr::mutate(
        w = rep(
            1,
            n()
        ),
        w = replace(
            w,
            5,
            Inf
        )
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_inf_w,
        weights = w
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_inf_w,
        weights = w
    ),
    error = function(e) e
)

econ_result
lm_result


# Negative weights -> as above
df_negative_w = mtcars %>%
    dplyr::mutate(
        w = rep(
            1,
            n()
        ),
        w = replace(
            w,
            5,
            -1
        )
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_negative_w,
        weights = w
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_negative_w,
        weights = w
    ),
    error = function(e) e
)

econ_result
lm_result


# NaN predictor -> warning in first part expected
df_nan = mtcars %>%
    dplyr::mutate(
        wt = replace(
            wt,
            5,
            NaN
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_nan
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_nan
)

nobs(m_econ)
nobs(m_lm)

m_econ$used
as.integer(row.names(model.frame(m_lm)))

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    vcov(m_econ),
    vcov(m_lm)
)

match(
    row.names(model.frame(m_lm)),
    row.names(mtcars)
)

m_econ$used

match(
    row.names(model.frame(m_lm)),
    row.names(mtcars)
)

identical(
    m_econ$used,
    match(
        row.names(model.frame(m_lm)),
        row.names(mtcars)
    )
)


# NaN weight -> warning edxpected
df_nan_w = mtcars %>%
    dplyr::mutate(
        w = rep(
            1,
            n()
        ),
        w = replace(
            w,
            5,
            NaN
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_nan_w,
    weights = w
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_nan_w,
    weights = w
)

nobs(m_econ)
nobs(m_lm)

m_econ$used

match(
    row.names(model.frame(m_lm)),
    row.names(mtcars)
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    vcov(m_econ),
    vcov(m_lm)
)


# NaN outcome
df_nan_y = mtcars %>%
    dplyr::mutate(
        mpg = replace(
            mpg,
            5,
            NaN
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_nan_y
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_nan_y
)

nobs(m_econ)
nobs(m_lm)

m_econ$used

match(
    row.names(model.frame(m_lm)),
    row.names(mtcars)
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    fitted(m_econ),
    fitted(m_lm)
)

all.equal(
    vcov(m_econ),
    vcov(m_lm)
)


# Extremely tiny positive weights (stress test) -> warnings wanted
df_tiny_w = mtcars %>%
    dplyr::mutate(
        w_tiny = rep(
            1e-300,
            n()
        )
    )

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_tiny_w,
    weights = w_tiny
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_tiny_w,
    weights = w_tiny
)

nobs(m_econ)
nobs(m_lm)

m_econ$rank
m_lm$rank

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-10
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-10
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-10
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-10
)


# Extremely large weights
df_huge_w = mtcars %>%
    dplyr::mutate(
        w_huge = rep(
            1e300,
            n()
        )
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_huge_w,
        weights = w_huge
    ),
    warning = function(w) w,
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_huge_w,
        weights = w_huge
    ),
    warning = function(w) w,
    error = function(e) e
)

econ_result
lm_result


# numerical comparison
m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_huge_w,
    weights = w_huge
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_huge_w,
    weights = w_huge
)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-10
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-10
)

all.equal(
    residuals(m_econ),
    residuals(m_lm),
    tolerance = 1e-10
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-10
)


# Extremely small nonzero predictor variation
df_precision = mtcars %>%
    dplyr::mutate(
        wt_precision = wt + seq_along(wt) * 1e-14
    )

m_econ = econ_lm(
    mpg ~ wt + wt_precision + hp,
    data = df_precision
)

m_lm = lm(
    mpg ~ wt + wt_precision + hp,
    data = df_precision
)

m_econ$rank
m_lm$rank

coef(m_econ)
coef(m_lm)

all.equal(
    coef(m_econ),
    coef(m_lm),
    tolerance = 1e-10
)

all.equal(
    fitted(m_econ),
    fitted(m_lm),
    tolerance = 1e-10
)

all.equal(
    vcov(m_econ),
    vcov(m_lm),
    tolerance = 1e-10
)


#### 27. API and Error Validation ####

# entirely omitted data -> plain text errors
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp
    ),
    error = function(e) e
)

econ_result
lm_result


# object data is not a data frame -> as above
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = 1:10
    ),
    error = function(e) e
)

econ_result


# Matrix as data
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = as.matrix(mtcars)
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = as.matrix(mtcars)
    ),
    error = function(e) e
)

econ_result
lm_result


# confirming tibble support before fixing
df_tibble = tibble::as_tibble(
    mtcars
)

is.data.frame(df_tibble)

m_econ = econ_lm(
    mpg ~ wt + hp,
    data = df_tibble
)

m_lm = lm(
    mpg ~ wt + hp,
    data = df_tibble
)

all.equal(
    coef(m_econ),
    coef(m_lm)
)

all.equal(
    vcov(m_econ),
    vcov(m_lm)
)


# Invalid vcov
tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        vcov = "banana"
    ),
    error = function(e) e
)


# Cluster covariance without cluster
tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        vcov = "CR2"
    ),
    error = function(e) e
)


# cluster supplied with non-clustered covariance
tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        vcov = "HC1",
        cluster = cyl
    ),
    error = function(e) e
)


# Invalid cluster length
bad_cluster = 1:10

tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        vcov = "CR2",
        cluster = bad_cluster
    ),
    error = function(e) e
)


# Invalid weights length
bad_weights = 1:10

tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        weights = bad_weights
    ),
    error = function(e) e
)


# Non-numeric weights
df_bad_weights = mtcars %>%
    dplyr::mutate(
        w = rep(
            "one",
            n()
        )
    )

tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_bad_weights,
        weights = w
    ),
    error = function(e) e
)


# Logical weights before changing validation
df_logical_w = mtcars %>%
    dplyr::mutate(
        w = rep(
            c(
                TRUE,
                FALSE
            ),
            length.out = n()
        )
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_logical_w,
        weights = w
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_logical_w,
        weights = w
    ),
    error = function(e) e
)

econ_result
lm_result


# factor weights
df_factor_w = mtcars %>%
    dplyr::mutate(
        w = factor(
            rep(
                c(
                    1,
                    2
                ),
                length.out = n()
            )
        )
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_factor_w,
        weights = w
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_factor_w,
        weights = w
    ),
    error = function(e) e
)

econ_result
lm_result


# all zero weights
df_zero_w = mtcars %>%
    dplyr::mutate(
        w = 0
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_zero_w,
        weights = w
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_zero_w,
        weights = w
    ),
    error = function(e) e
)

econ_result
lm_result


# empty df
df_empty = mtcars %>%
    dplyr::slice(
        0
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_empty
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_empty
    ),
    error = function(e) e
)

econ_result
lm_result


# All observations removed by missingness
df_all_missing = mtcars %>%
    dplyr::mutate(
        wt = NA_real_
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_all_missing
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_all_missing
    ),
    error = function(e) e
)

econ_result
lm_result


# Missing outcome variable
econ_result = tryCatch(
    econ_lm(
        nonexistent_y ~ wt + hp,
        data = mtcars
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        nonexistent_y ~ wt + hp,
        data = mtcars
    ),
    error = function(e) e
)

econ_result
lm_result


# Missing predictor variable
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + nonexistent_x,
        data = mtcars
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + nonexistent_x,
        data = mtcars
    ),
    error = function(e) e
)

econ_result
lm_result


# Invalid formula object
econ_result = tryCatch(
    econ_lm(
        "mpg ~ wt + hp",
        data = mtcars
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        "mpg ~ wt + hp",
        data = mtcars
    ),
    error = function(e) e
)

econ_result
lm_result


# Multi-response formula
econ_result = tryCatch(
    econ_lm(
        cbind(mpg, disp) ~ wt + hp,
        data = mtcars
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        cbind(mpg, disp) ~ wt + hp,
        data = mtcars
    ),
    error = function(e) e
)

econ_result
lm_result


# Non-numeric response -> plain text warning
df_factor_y = mtcars %>%
    dplyr::mutate(
        mpg_factor = factor(
            mpg > median(mpg)
        )
    )

econ_result = tryCatch(
    econ_lm(
        mpg_factor ~ wt + hp,
        data = df_factor_y
    ),
    warning = function(w) w,
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg_factor ~ wt + hp,
        data = df_factor_y
    ),
    warning = function(w) w,
    error = function(e) e
)

econ_result
lm_result


# character response
df_character_y = mtcars %>%
    dplyr::mutate(
        mpg_character = as.character(
            mpg
        )
    )

tryCatch(
    econ_lm(
        mpg_character ~ wt + hp,
        data = df_character_y
    ),
    error = function(e) e
)


# One-sided formula / no response
econ_result = tryCatch(
    econ_lm(
        ~ wt + hp,
        data = mtcars
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        ~ wt + hp,
        data = mtcars
    ),
    error = function(e) e
)

econ_result
lm_result


# Invalid subset length
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = c(TRUE, FALSE)
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = c(TRUE, FALSE)
    ),
    error = function(e) e
)

econ_result
lm_result


# Numeric subset
econ_result = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = c(1, 3, 5, 10, 20)
)

lm_result = lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = c(1, 3, 5, 10, 20)
)

nobs(econ_result)
nobs(lm_result)

all.equal(
    coef(econ_result),
    coef(lm_result)
)

all.equal(
    fitted(econ_result),
    fitted(lm_result)
)

all.equal(
    vcov(econ_result),
    vcov(lm_result)
)


# Subset containing NA
subset_na = rep(
    TRUE,
    nrow(mtcars)
)

subset_na[c(5, 10)] = NA

econ_result = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = subset_na
)

lm_result = lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = subset_na
)

nobs(econ_result)
nobs(lm_result)

all.equal(
    coef(econ_result),
    coef(lm_result)
)

all.equal(
    fitted(econ_result),
    fitted(lm_result)
)

all.equal(
    vcov(econ_result),
    vcov(lm_result)
)


# empty subset
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = FALSE
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = FALSE
    ),
    error = function(e) e
)

econ_result
lm_result


# Invalid na.action
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        na.action = "banana"
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = mtcars,
        na.action = "banana"
    ),
    error = function(e) e
)

econ_result
lm_result


# na.action supplied as a character name
df_na_action = mtcars %>%
    dplyr::mutate(
        hp = replace(
            hp,
            5,
            NA
        )
    )

econ_result = econ_lm(
    mpg ~ wt + hp,
    data = df_na_action,
    na.action = "na.omit"
)

lm_result = lm(
    mpg ~ wt + hp,
    data = df_na_action,
    na.action = "na.omit"
)

nobs(econ_result)
nobs(lm_result)

all.equal(
    coef(econ_result),
    coef(lm_result)
)

all.equal(
    fitted(econ_result),
    fitted(lm_result)
)

all.equal(
    vcov(econ_result),
    vcov(lm_result)
)


# Custom na.action function
econ_result = econ_lm(
    mpg ~ wt + hp,
    data = df_na_action,
    na.action = stats::na.omit
)

lm_result = lm(
    mpg ~ wt + hp,
    data = df_na_action,
    na.action = stats::na.omit
)

nobs(econ_result)
nobs(lm_result)

all.equal(
    coef(econ_result),
    coef(lm_result)
)

all.equal(
    fitted(econ_result),
    fitted(lm_result)
)

all.equal(
    vcov(econ_result),
    vcov(lm_result)
)


# single cluster
df_one_cluster = mtcars %>%
    dplyr::mutate(
        cluster_id = 1
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_one_cluster,
        vcov = "CR2",
        cluster = cluster_id
    ),
    error = function(e) e,
    warning = function(w) w
)

econ_result


# All cluster IDs missing
df_missing_cluster = mtcars %>%
    dplyr::mutate(
        cluster_id = NA_integer_
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_missing_cluster,
        vcov = "CR2",
        cluster = cluster_id
    ),
    error = function(e) e
)

econ_result


# Non-vector cluster input
df_cluster_matrix = mtcars %>%
    dplyr::mutate(
        cluster_a = rep(1:8, each = 4),
        cluster_b = rep(1:4, each = 8)
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_cluster_matrix,
        vcov = "CR2",
        cluster = cbind(
            cluster_a,
            cluster_b
        )
    ),
    error = function(e) e,
    warning = function(w) w
)

econ_result


# Missing cluster variable
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        vcov = "CR2",
        cluster = nonexistent_cluster
    ),
    error = function(e) e
)

econ_result


# Missing weights variable
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        weights = nonexistent_weights
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = mtcars,
        weights = nonexistent_weights
    ),
    error = function(e) e
)

econ_result
lm_result


# Invalid subset variable
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = nonexistent_subset
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = nonexistent_subset
    ),
    error = function(e) e
)

econ_result
lm_result


# Invalid subset type
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = "banana"
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = "banana"
    ),
    error = function(e) e
)

econ_result
lm_result


# Complex weights
df_complex_weights = mtcars %>%
    dplyr::mutate(
        w_complex = rep(
            1 + 1i,
            n()
        )
    )

econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = df_complex_weights,
        weights = w_complex
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = df_complex_weights,
        weights = w_complex
    ),
    error = function(e) e
)

econ_result
lm_result


# NA in numeric subset indices
econ_result = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = c(1, 3, NA, 5, 10)
)

lm_result = lm(
    mpg ~ wt + hp,
    data = mtcars,
    subset = c(1, 3, NA, 5, 10)
)

nobs(econ_result)
nobs(lm_result)

all.equal(
    coef(econ_result),
    coef(lm_result)
)

all.equal(
    fitted(econ_result),
    fitted(lm_result)
)

all.equal(
    vcov(econ_result),
    vcov(lm_result)
)


# Zero-length subset
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = logical(0)
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = logical(0)
    ),
    error = function(e) e
)

econ_result
lm_result


# Zero-length weights
econ_result = tryCatch(
    econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        weights = numeric(0)
    ),
    error = function(e) e
)

lm_result = tryCatch(
    lm(
        mpg ~ wt + hp,
        data = mtcars,
        weights = numeric(0)
    ),
    error = function(e) e
)

econ_result
lm_result


#### 28. compatibility and final regression testing ####

# broom::tidy() compatibility
a = econ_lm(
    mpg ~ wt + hp,
    data = mtcars
)

b = lm(
    mpg ~ wt + hp,
    data = mtcars
)

broom::tidy(a)
broom::tidy(b)


# broom::tidy(..., conf.int = FALSE)
broom::tidy(
    a,
    conf.int = FALSE
)


# custom confidence level
broom::tidy(
    a,
    conf.int = TRUE,
    conf.level = 0.90
)


# broom::glance() compatibility
broom::glance(a)
broom::glance(b)


# AIC/BIC direct comparison
c(
    econ_aic = broom::glance(a)$aic,
    lm_aic = AIC(b),
    difference = AIC(b) - broom::glance(a)$aic
)

c(
    econ_bic = broom::glance(a)$bic,
    lm_bic = BIC(b),
    difference = BIC(b) - broom::glance(a)$bic
)


# modelsummary compatibility
modelsummary::modelsummary(
    list(
        "EconR" = a,
        "Base R" = b
    ),
    output = "data.frame"
)


# sandwich external compatibility
a_hc3 = sandwich::vcovHC(
    a,
    type = "HC3"
)

b_hc3 = sandwich::vcovHC(
    b,
    type = "HC3"
)

all.equal(
    unname(a_hc3),
    unname(b_hc3),
    tolerance = 1e-10
)

sandwich::bread(a)
sandwich::estfun(a)[1:5, ]


# update() compatibility
a_updated = update(
    a,
    . ~ . + qsec
)

b_updated = update(
    b,
    . ~ . + qsec
)

coef(a_updated)
coef(b_updated)

all.equal(
    coef(a_updated),
    coef(b_updated),
    tolerance = 1e-10
)


# formula(), terms(), and model.frame()
formula(a)
formula(b)

terms(a)
terms(b)

model.frame(a)[1:5, ]
model.frame(b)[1:5, ]


# model.matrix() and weights()
all.equal(
    model.matrix(a),
    model.matrix(b),
    tolerance = 1e-10
)

head(weights(a))
head(weights(b))

nobs(a)
nobs(b)


# predict() compatibility through generic dispatch
all.equal(
    predict(a),
    predict(b),
    tolerance = 1e-10
)

newdata = data.frame(
    wt = c(2, 3, 4),
    hp = c(80, 150, 220)
)

all.equal(
    predict(a, newdata = newdata),
    predict(b, newdata = newdata),
    tolerance = 1e-10
)


# residuals(), fitted(), coef(), vcov()
all.equal(
    residuals(a),
    residuals(b),
    tolerance = 1e-10
)

all.equal(
    fitted(a),
    fitted(b),
    tolerance = 1e-10
)

all.equal(
    coef(a),
    coef(b),
    tolerance = 1e-10
)

all.equal(
    vcov(a),
    vcov(b),
    tolerance = 1e-10
)


# na.exclude through downstream generics
d = mtcars %>%
    dplyr::mutate(
        hp = replace(
            hp,
            c(5, 10),
            NA
        )
    )

a_ex = econ_lm(
    mpg ~ wt + hp,
    data = d,
    na.action = stats::na.exclude
)

b_ex = lm(
    mpg ~ wt + hp,
    data = d,
    na.action = stats::na.exclude
)

length(residuals(a_ex))
length(residuals(b_ex))

which(is.na(residuals(a_ex)))
which(is.na(residuals(b_ex)))

all.equal(
    residuals(a_ex),
    residuals(b_ex),
    tolerance = 1e-10
)


# weighted-model compatibility
d_w = mtcars %>%
    dplyr::mutate(
        w = seq(0.5, 2, length.out = n())
    )

a_w = econ_lm(
    mpg ~ wt + hp,
    data = d_w,
    weights = w
)

b_w = lm(
    mpg ~ wt + hp,
    data = d_w,
    weights = w
)

all.equal(
    coef(a_w),
    coef(b_w),
    tolerance = 1e-10
)

all.equal(
    vcov(a_w),
    vcov(b_w),
    tolerance = 1e-10
)

all.equal(
    fitted(a_w),
    fitted(b_w),
    tolerance = 1e-10
)

nobs(a_w)
nobs(b_w)


# rank-deficient model through standard generics
d_rank = mtcars %>%
    dplyr::mutate(
        wt_copy = wt
    )

a_rank = econ_lm(
    mpg ~ wt + wt_copy + hp,
    data = d_rank
)

b_rank = lm(
    mpg ~ wt + wt_copy + hp,
    data = d_rank
)

coef(a_rank)
coef(b_rank)

all.equal(
    coef(a_rank),
    coef(b_rank),
    tolerance = 1e-10
)

all.equal(
    fitted(a_rank),
    fitted(b_rank),
    tolerance = 1e-10
)


# rank-deficient vcov() and summary()
vcov(a_rank)
vcov(b_rank)

summary(a_rank)$coefficients
summary(b_rank)$coefficients


# broom::tidy() with rank deficiency
broom::tidy(a_rank)
broom::tidy(b_rank)


# modelsummary with rank deficiency
modelsummary::modelsummary(
    list(
        "EconR" = a_rank,
        "Base R" = b_rank
    ),
    output = "data.frame"
)


# confint() compatibility
confint(a)
confint(b)

# first diagnose the dispatch
getS3method(
    "confint",
    "econ_lm",
    optional = TRUE
)

getS3method(
    "confint",
    "default"
)


# first check lm() parameter-selection behavior
confint(
    b,
    parm = c("wt", "hp"),
    level = 0.90
)

confint(
    b,
    parm = 2,
    level = 0.90
)


# named parameter selection
confint(
    a,
    parm = c("wt", "hp"),
    level = 0.90
)

confint(
    b,
    parm = c("wt", "hp"),
    level = 0.90
)


# numeric parameter selection
confint(
    a,
    parm = 2,
    level = 0.90
)

confint(
    b,
    parm = 2,
    level = 0.90
)


# rank-deficient confidence intervals
confint(a_rank)
confint(b_rank)


# weighted confidence intervals
confint(a_w)
confint(b_w)

all.equal(
    confint(a_w),
    confint(b_w),
    tolerance = 1e-10
)


# anova() compatibility
anova(a)
anova(b)


# wls anova
anova(a_w)
anova(b_w)


# factor anova
a_factor = econ_lm(
    mpg ~ factor(cyl) + wt,
    data = mtcars
)

anova(a_factor)
anova(b_factor)


# interaction anova
a_interaction = econ_lm(
    mpg ~ factor(cyl) * wt,
    data = mtcars
)

b_interaction = lm(
    mpg ~ factor(cyl) * wt,
    data = mtcars
)

anova(a_interaction)
anova(b_interaction)


# rank deficiency anova
anova(a_rank)
anova(b_rank)


# offset anova
a_off = econ_lm(
    mpg ~ wt + hp + offset(qsec),
    data = mtcars
)

anova(a_off)
anova(b_off)


# no-intercept ANOVA
a_no_intercept = econ_lm(
    mpg ~ 0 + wt + hp,
    data = mtcars
)

b_no_intercept = lm(
    mpg ~ 0 + wt + hp,
    data = mtcars
)

anova(a_no_intercept)
anova(b_no_intercept)


# anova w zero weights
w_zero = rep(
    1,
    nrow(mtcars)
)

w_zero[c(2, 7, 15)] = 0

a_zero = econ_lm(
    mpg ~ wt + hp,
    data = mtcars,
    weights = w_zero
)

b_zero = lm(
    mpg ~ wt + hp,
    data = mtcars,
    weights = w_zero
)

anova(a_zero)
anova(b_zero)


# anova w na.exclude
mt_na = mtcars %>%
    mutate(
        mpg = replace(
            mpg,
            c(5, 10),
            NA
        )
    )

a_na_anova = econ_lm(
    mpg ~ wt + hp,
    data = mt_na,
    na.action = na.exclude
)

b_na_anova = lm(
    mpg ~ wt + hp,
    data = mt_na,
    na.action = na.exclude
)

anova(a_na_anova)
anova(b_na_anova)


# inspect exact multiple model structure
str(
    anova(
        b_small,
        b
    )
)


# model w diff estimation samples
mt_diff = mtcars %>%
    mutate(
        hp = replace(
            hp,
            1,
            NA
        )
    )

b_diff_1 = lm(
    mpg ~ wt,
    data = mt_diff
)

b_diff_2 = lm(
    mpg ~ wt + hp,
    data = mt_diff
)

anova(
    b_diff_1,
    b_diff_2
)


# same sample size diff obs
mt_diff_1 = mtcars %>%
    mutate(
        mpg = replace(
            mpg,
            1,
            NA
        )
    )

mt_diff_2 = mtcars %>%
    mutate(
        mpg = replace(
            mpg,
            2,
            NA
        )
    )

b_diff_1 = lm(
    mpg ~ wt,
    data = mt_diff_1
)

b_diff_2 = lm(
    mpg ~ wt + hp,
    data = mt_diff_2
)

anova(
    b_diff_1,
    b_diff_2
)


# reversing model order
anova(
    a,
    a_small
)

anova(
    b,
    b_small
)


# 3 model comparisson
a_0 = econ_lm(
    mpg ~ 1,
    data = mtcars
)

anova(
    a_0,
    a_small,
    a
)


# one last source check
getAnywhere("anova.lmlist")


# basic nested model comparisson
a_small = econ_lm(
    mpg ~ wt,
    data = mtcars
)

anova(
    a_small,
    a
)


# different estimation sample sizes
a_diff_1 = econ_lm(
    mpg ~ wt,
    data = mt_diff
)

a_diff_2 = econ_lm(
    mpg ~ wt + hp,
    data = mt_diff
)

anova(
    a_diff_1,
    a_diff_2
)


# same size, different observations
a_diff_same_1 = econ_lm(
    mpg ~ wt,
    data = mt_diff_1
)

a_diff_same_2 = econ_lm(
    mpg ~ wt + hp,
    data = mt_diff_2
)

anova(
    a_diff_same_1,
    a_diff_same_2
)


# diff response variables
a_other_response = econ_lm(
    qsec ~ wt + hp,
    data = mtcars
)

anova(
    a,
    a_other_response
)


# explicit scale
anova(
    a_small,
    a,
    scale = 10
)

anova(
    b_small,
    b,
    scale = 10
)


# test = NULL
anova(
    a_small,
    a,
    test = NULL
)

anova(
    b_small,
    b,
    test = NULL
)








































#### Dev Tool Testing ####
devtools::load_all()
devtools::test()

devtools::document()

devtools::check()

