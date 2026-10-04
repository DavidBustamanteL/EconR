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

