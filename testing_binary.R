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


###### 2. Logit Fitted Probabilities ####

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












































#### Dev Tool Testing ####
devtools::load_all()
devtools::test()

devtools::document()

devtools::check()