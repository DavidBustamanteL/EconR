################################# Testing Multinomial Models EconR ###############################
##################################################################################################

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
library(nnet)
library(modelsummary)

# Installing test pack if required
if (!requireNamespace("EconR", quietly = TRUE)) {
    install.packages(
        "G:/My Drive/5. CURSOS ONLINE/R Practice FULL/Econ Own Package",
        repos = NULL,
        type = "source"
    )
}

library(EconR)


###### 1. Verify EconR Summary Output ####

s = summary(model)

print(s)

EconR:::print.summary.econ_multinomial(s)


###### 2. Verify HC1 Standard Errors ####

s$standard.errors

matrix(
    sqrt(diag(vcov(model))),
    nrow = nrow(coef(model)),
    byrow = TRUE,
    dimnames = dimnames(coef(model))
)

all.equal(
    s$standard.errors,
    matrix(
        sqrt(diag(vcov(model))),
        nrow = nrow(coef(model)),
        byrow = TRUE,
        dimnames = dimnames(coef(model))
    )
)




###### 3. Reproduce the Failing Model ####

dat = make_multinomial_sample(n = 1000, seed = 1238)

ref = nnet::multinom(
    y ~ x1 + x2,
    data = dat,
    weights = w,
    subset = keep,
    Hess = TRUE,
    trace = FALSE
)

fit = econ_multinomial(
    y ~ x1 + x2,
    data = dat,
    weights = w,
    subset = keep,
    vcov = "CR0",
    cluster = cluster_id
)

c(
    n_ref = nrow(model.frame(ref)),
    n_econ = nobs(fit),
    n_cluster = length(fit$cluster),
    n_weights = length(fit$econ_weights)
)

all.equal(unname(coef(ref)), unname(coef(fit)))


###### 4. Compare Weighted Scores ####

S_ref = sandwich::estfun(ref)
S_econ = sandwich::estfun(fit)

dim(S_ref)
dim(S_econ)

all.equal(
    unname(S_ref),
    unname(S_econ),
    tolerance = 1e-7
)


###### 5. Compare Covariance Components ####

V_ref = vcov(ref)

V_econ_classical = getS3method("vcov", "multinom")(fit)

all.equal(
    unname(V_ref),
    unname(V_econ_classical),
    tolerance = 1e-8
)

U = rowsum(
    S_ref,
    group = fit$cluster,
    reorder = TRUE
)

V_manual = V_ref %*% crossprod(U) %*% V_ref

V_sandwich = sandwich::vcovCL(
    ref,
    cluster = fit$cluster,
    type = "HC0",
    cadjust = FALSE
)

all.equal(
    unname(vcov(fit)),
    unname(V_manual),
    tolerance = 1e-7
)

all.equal(
    unname(V_manual),
    unname(V_sandwich),
    tolerance = 1e-7
)

cbind(
    EconR = diag(vcov(fit)),
    Manual = diag(V_manual),
    Sandwich = diag(V_sandwich)
)

























#### Dev Tool Testing ####
devtools::load_all()
devtools::test(filter = "multinomial")

devtools::document()

devtools::check()
