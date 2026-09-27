################################# Further Testing Package ##################################
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
library(broom)
library(estimatr)
library(WeightIt)
library(sandwich)
library(clubSandwich)
library(modelsummary)

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

# WeightIt -> HC0 in WeightIT is HC1 in all other packages
model_wei = lm_weightit(wage ~ educ + exper + nonwhite + female + married, data = wages, vcov = "HC0")
summary(model_wei)
tidy(model_wei)





# Interactions
model1 = lm_robust(
    wage ~ educ + exper * female,
    data = wages, se_type = "HC1")
tidy(model1)

model_econr1 = econ_lm(
    wage ~ educ + exper * female,
    data = wages, vcov = "HC1")
tidy(model_econr1)
