################################# Further Testing Package ##################################
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
library(broom)
library(estimatr)
library(WeightIt)
library(sandwich)
library(clubSandwich)
library(modelsummary)
library(lmtest)

library(EconR)

# importing the data
wages = read.csv("wage1.csv", header = TRUE, sep = ",")


#### 1. Standard OLS ####

# estimatr
model = lm_robust(wage ~ educ + exper + nonwhite + female + married, data = wages, se_type = "HC1")
summary(model)
tidy(model)

# EconR
model_econr = econ_lm(wage ~ educ + exper + nonwhite + female + married, data = wages, vcov = "HC1")
summary(model_econr)
tidy(model_econr)
modelsummary(
    model_econr,
    stars = TRUE,
    output = "default"
)

# WeightIt -> HC0 in WeightIT is HC1 in all other packages
model_wei = lm_weightit(wage ~ educ + exper + nonwhite + female + married, data = wages, vcov = "HC0")
summary(model_wei)
tidy(model_wei)

# model summary expanded
model_1 = econ_lm(
    wage ~ educ,
    data = wages,
    vcov = "HC1"
)

model_2 = econ_lm(
    wage ~ educ + exper,
    data = wages,
    vcov = "HC1"
)

model_3 = econ_lm(
    wage ~ educ + exper + nonwhite + female + married,
    data = wages,
    vcov = "HC1"
)

model_4 = econ_lm(
    wage ~ educ * female + exper + nonwhite + married,
    data = wages,
    vcov = "HC1"
)

models = list(
    "Education" = model_1,
    "Human Capital" = model_2,
    "Full Model" = model_3,
    "Interaction" = model_4
)

modelsummary(
    models,
    stars = TRUE,
    output = "default"
)


# Interactions
model1 = lm_robust(
    wage ~ educ + exper * female,
    data = wages, se_type = "HC1")
tidy(model1)

model_econr1 = econ_lm(
    wage ~ educ + exper * female,
    data = wages, vcov = "HC1")
tidy(model_econr1)

car::linearHypothesis(
    model_econr1,
    c("female = 0", "exper:female = 0"),
    vcov. = vcov(model_econr1)
)
lmtest::waldtest(model_econr1, c("female", "exper:female"))
waldtest(model_econr1, . ~ . - exper - exper:female, test = "F")


# binary #

#### 1. Standard Logit ####

# glm
model = glm(
    married ~ educ + exper + nonwhite + female + wage,
    data = wages,
    family = binomial(link = "logit")
)

summary(model)
tidy(model)


# EconR
model_econr = econ_binary(
    married ~ educ + exper + nonwhite + female + wage,
    data = wages,
    link = "logit",
    vcov = "classical"
)

summary(model_econr)
tidy(model_econr)


#### 2. Compare Results ####

# Coefficients
all.equal(
    coef(model),
    coef(model_econr),
    tolerance = 1e-10
)

# Covariance matrices
all.equal(
    vcov(model),
    vcov(model_econr),
    tolerance = 1e-10
)

# Predicted probabilities
all.equal(
    predict(model, type = "response"),
    predict(model_econr, type = "response"),
    tolerance = 1e-10
)


#### 3. Robust Standard Errors (HC1) ####

model_econr_HC1 = econ_binary(
    married ~ educ + exper + nonwhite + female + wage,
    data = wages,
    link = "logit",
    vcov = "HC1"
)

# Compare standard errors
cbind(
    glm = sqrt(diag(sandwich::vcovHC(model, type = "HC1"))),
    EconR = sqrt(diag(vcov(model_econr_HC1)))
)


#### 4. Modelsummary ####

modelsummary(
    list(
        "glm" = model,
        "EconR" = model_econr_HC1
    ),
    vcov = list(
        sandwich::vcovHC(model, type = "HC1"),
        vcov(model_econr_HC1)
    ),
    stars = TRUE,
    output = "default"
)



#### 5. Logit with Interaction ####

# glm
model_int = glm(
    married ~ educ * female + exper + nonwhite + wage,
    data = wages,
    family = binomial(link = "logit")
)

summary(model_int)
tidy(model_int)


# EconR
model_econr_int = econ_binary(
    married ~ educ * female + exper + nonwhite + wage,
    data = wages,
    link = "logit",
    vcov = "HC1"
)

summary(model_econr_int)
tidy(model_econr_int)


# Compare coefficients
all.equal(
    coef(model_int),
    coef(model_econr_int),
    tolerance = 1e-10
)

# Compare HC1 covariance matrices
all.equal(
    sandwich::vcovHC(model_int, type = "HC1"),
    vcov(model_econr_int),
    tolerance = 1e-10
)

# Compare predicted probabilities
all.equal(
    predict(model_int, type = "response"),
    predict(model_econr_int, type = "response"),
    tolerance = 1e-10
)


#### 6. Logit with Quadratic Term ####

# glm
model_quad = glm(
    married ~ educ + exper + I(exper^2) + nonwhite + female + wage,
    data = wages,
    family = binomial(link = "logit")
)

summary(model_quad)
tidy(model_quad)


# EconR
model_econr_quad = econ_binary(
    married ~ educ + exper + I(exper^2) + nonwhite + female + wage,
    data = wages,
    link = "logit",
    vcov = "HC1"
)

summary(model_econr_quad)
tidy(model_econr_quad)


# Compare coefficients
all.equal(
    coef(model_quad),
    coef(model_econr_quad),
    tolerance = 1e-10
)

# Compare HC1 covariance matrices
all.equal(
    sandwich::vcovHC(model_quad, type = "HC1"),
    vcov(model_econr_quad),
    tolerance = 1e-10
)

# Compare predicted probabilities
all.equal(
    predict(model_quad, type = "response"),
    predict(model_econr_quad, type = "response"),
    tolerance = 1e-10
)


#### 7. Combined Interaction and Quadratic Terms ####

# glm
model_full = glm(
    married ~ educ * female + exper + I(exper^2) + nonwhite + wage,
    data = wages,
    family = binomial(link = "logit")
)

summary(model_full)
tidy(model_full)


# EconR
model_econr_full = econ_binary(
    married ~ educ * female + exper + I(exper^2) + nonwhite + wage,
    data = wages,
    link = "logit",
    vcov = "HC1"
)

summary(model_econr_full)
tidy(model_econr_full)


# Compare coefficients
all.equal(
    coef(model_full),
    coef(model_econr_full),
    tolerance = 1e-10
)

# Compare HC1 covariance matrices
all.equal(
    sandwich::vcovHC(model_full, type = "HC1"),
    vcov(model_econr_full),
    tolerance = 1e-10
)

# Compare predicted probabilities
all.equal(
    predict(model_full, type = "response"),
    predict(model_econr_full, type = "response"),
    tolerance = 1e-10
)


#### 8. Modelsummary ####

modelsummary(
    list(
        "Interaction (glm)" = model_int,
        "Interaction (EconR)" = model_econr_int,
        "Quadratic (glm)" = model_quad,
        "Quadratic (EconR)" = model_econr_quad,
        "Combined (glm)" = model_full,
        "Combined (EconR)" = model_econr_full
    ),
    vcov = list(
        sandwich::vcovHC(model_int, type = "HC1"),
        vcov(model_econr_int),
        sandwich::vcovHC(model_quad, type = "HC1"),
        vcov(model_econr_quad),
        sandwich::vcovHC(model_full, type = "HC1"),
        vcov(model_econr_full)
    ),
    stars = TRUE,
    output = "default"
)


