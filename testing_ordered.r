############################### Testing Ordered Models EconR ###############################
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
library(MASS)
library(rms)
library(ordinal)

# Installing test pack if required
if (!requireNamespace("EconR", quietly = TRUE)) {
    install.packages(
        "G:/My Drive/5. CURSOS ONLINE/R Practice FULL/Econ Own Package",
        repos = NULL,
        type = "source"
    )
}

library(EconR)


#### 1. Ordered Model Reference Validation ####

set.seed(123)

n = 1000

df = data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
)

# Latent outcome
eta = 0.8 * df$x1 - 0.5 * df$x2 + 0.3 * df$x3

# Ordered outcome: Low < Medium < High
latent = eta + rlogis(n)

df$y = cut(
    latent,
    breaks = c(-Inf, -0.5, 0.7, Inf),
    labels = c("Low", "Medium", "High"),
    ordered_result = TRUE
)

# Ordered Logit
ref_logit = MASS::polr(
    y ~ x1 + x2 + x3,
    data = df,
    method = "logistic",
    Hess = TRUE
)

# Ordered Probit
ref_probit = MASS::polr(
    y ~ x1 + x2 + x3,
    data = df,
    method = "probit",
    Hess = TRUE
)

# Inspect results
summary(ref_logit)
summary(ref_probit)

# Thresholds
ref_logit$zeta
ref_probit$zeta

# Predicted category probabilities
head(predict(ref_logit, type = "probs"))
head(predict(ref_probit, type = "probs"))


#### 2. Validate Ordered Estimation ####

econ_logit = econ_ordered(
    y ~ x1 + x2 + x3,
    data = df,
    link = "logit"
)

econ_probit = econ_ordered(
    y ~ x1 + x2 + x3,
    data = df,
    link = "probit"
)

# Compare slope coefficients
all.equal(
    coef(econ_logit),
    coef(ref_logit)
)

all.equal(
    coef(econ_probit),
    coef(ref_probit)
)

# Compare thresholds
all.equal(
    econ_logit$zeta,
    ref_logit$zeta
)

all.equal(
    econ_probit$zeta,
    ref_probit$zeta
)

# Compare predicted probabilities
all.equal(
    predict(econ_logit, type = "probs"),
    predict(ref_logit, type = "probs")
)

all.equal(
    predict(econ_probit, type = "probs"),
    predict(ref_probit, type = "probs")
)


#### 3. Ordered Model Missing-Data Validation ####

set.seed(456)

df_na = df

# Introduce missing values
df_na$x1[sample(nrow(df_na), 50)] = NA
df_na$x2[sample(nrow(df_na), 70)] = NA
df_na$x3[sample(nrow(df_na), 40)] = NA
df_na$y[sample(nrow(df_na), 30)] = NA

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_na,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_na,
        link = link
    )

    cat("\nModel:", link, "\n")

    cat("Reference observations:", nobs(ref), "\n")
    cat("EconR observations:", nobs(econ), "\n")

    cat("Coefficients equal:",
        isTRUE(all.equal(coef(econ), coef(ref))), "\n")

    cat("Thresholds equal:",
        isTRUE(all.equal(econ$zeta, ref$zeta)), "\n")

    cat("Probabilities equal:",
        isTRUE(all.equal(
            predict(econ, type = "probs"),
            predict(ref, type = "probs")
        )), "\n")
}


#### 4. Ordered Models: Weights and Subsets ####

set.seed(789)

df_test = df

# Integer weights avoid fractional binomial-weight warnings
df_test$w = sample(1:3, nrow(df_test), replace = TRUE)
df_test$w[sample(nrow(df_test), 50)] = NA

# Subset indicator
df_test$keep = sample(c(TRUE, FALSE), nrow(df_test),
                    replace = TRUE, prob = c(0.8, 0.2))

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link
    )

    cat("\nModel:", link, "\n")

    cat("Reference observations:", nobs(ref), "\n")
    cat("EconR observations:", nobs(econ), "\n")

    cat("Coefficients equal:",
        isTRUE(all.equal(coef(econ), coef(ref))), "\n")

    cat("Thresholds equal:",
        isTRUE(all.equal(econ$zeta, ref$zeta)), "\n")

    cat("Prior weights equal:",
        isTRUE(all.equal(econ$weights, ref$weights)), "\n")
}


#### 5. Ordered Model Covariance Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link
    )

    V_ref = vcov(ref)
    V_econ = vcov(econ)

    cat("\nModel:", link, "\n")

    cat("Covariance matrices equal:",
        isTRUE(all.equal(V_econ, V_ref)), "\n")

    cat("Standard errors equal:",
        isTRUE(all.equal(
            sqrt(diag(V_econ)),
            sqrt(diag(V_ref))
        )), "\n")

    cat("Parameters in covariance matrix:",
        paste(rownames(V_econ), collapse = ", "), "\n")

    cat("Actual estimation rows:",
        length(econ$used), "\n")

    cat("Sum of prior weights:",
        sum(econ$weights), "\n")
}


#### 6. Ordered Model Weight Storage ####

devtools::load_all()

for (link in c("logit", "probit")) {

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link
    )

    cat("\nModel:", link, "\n")

    cat("Estimation rows:", length(econ$used), "\n")
    cat("Stored weights:", length(econ$weights), "\n")
    cat("Sum of weights:", sum(econ$weights), "\n")

    cat("Weights aligned:",
        length(econ$weights) == length(econ$used), "\n")

    cat("Covariance still valid:",
        isTRUE(all.equal(
            vcov(econ),
            vcov(MASS::polr(
                y ~ x1 + x2 + x3,
                data = df_test,
                weights = w,
                subset = keep,
                method = if (link == "logit") "logistic" else "probit",
                Hess = TRUE
            ))
        )), "\n")
}


#### 7. Ordered Model Score Validation ####

for (link in c("logit", "probit")) {

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link
    )

    S = sandwich::estfun(econ)

    cat("\nModel:", link, "\n")

    cat("Score dimensions:",
        paste(dim(S), collapse = " x "), "\n")

    cat("Score column names:",
        paste(colnames(S), collapse = ", "), "\n")

    cat("Scores finite:",
        all(is.finite(S)), "\n")

    cat("Score sums approximately zero:",
        max(abs(colSums(S))) < 1e-4, "\n")

    cat("Maximum absolute score sum:",
        max(abs(colSums(S))), "\n")
}


#### 8. Ordered Model Numerical Score Validation ####

for (link in c("logit", "probit")) {

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link
    )

    X = model.matrix(econ)
    X = X[, colnames(X) != "(Intercept)", drop = FALSE]

    y = as.integer(model.response(econ$model))
    weights = econ$weights

    K = length(coef(econ))
    theta = c(coef(econ), econ$zeta)

    # Weighted log-likelihood
    loglik = function(par) {

        beta = par[seq_len(K)]
        zeta = par[-seq_len(K)]

        eta = drop(X %*% beta)

        F = if (link == "logit") plogis else pnorm

        bounds = c(-Inf, zeta, Inf)

        p = F(bounds[y + 1L] - eta) -
            F(bounds[y] - eta)

        sum(weights * log(p))
    }

    # Central finite-difference gradient
    h = 1e-5

    numerical_gradient = sapply(
        seq_along(theta),
        function(j) {

            plus = theta
            minus = theta

            plus[j] = plus[j] + h
            minus[j] = minus[j] - h

            (loglik(plus) - loglik(minus)) / (2 * h)
        }
    )

    analytical_gradient = colSums(
        sandwich::estfun(econ)
    )

    cat("\nModel:", link, "\n")

    print(
        data.frame(
            parameter = names(theta),
            analytical = unname(analytical_gradient),
            numerical = unname(numerical_gradient),
            difference = unname(
                analytical_gradient - numerical_gradient
            )
        ),
        row.names = FALSE
    )

    cat(
        "Maximum absolute difference:",
        max(abs(
            analytical_gradient - numerical_gradient
        )),
        "\n"
    )
}


#### 9. Ordered Model HC0 Validation ####

for (link in c("logit", "probit")) {

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link,
        vcov = "HC0"
    )

    S = sandwich::estfun(econ)
    H = econ$Hessian

    V_manual = solve(H) %*%
        crossprod(S) %*%
        solve(H)

    V_econ = vcov(econ)

    cat("\nModel:", link, "\n")

    cat("HC0 matches manual calculation:",
        isTRUE(all.equal(
            unname(V_econ),
            unname(V_manual)
        )), "\n")

    cat("Covariance symmetric:",
        isTRUE(all.equal(
            unname(V_econ),
            unname(t(V_econ))
        )), "\n")

    cat("Standard errors finite:",
        all(is.finite(sqrt(diag(V_econ)))), "\n")

    cat("HC0 standard errors:\n")
    print(sqrt(diag(V_econ)))
}


#### 10. Ordered Model Independent HC0 Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link,
        vcov = "HC0"
    )

    V_econ = vcov(econ)

    V_ref = tryCatch(
        sandwich::sandwich(ref),
        error = function(e) e
    )

    cat("\nModel:", link, "\n")

    if (inherits(V_ref, "error")) {

        cat("Reference unavailable:",
            conditionMessage(V_ref), "\n")

    } else {

        cat("Independent covariance matches:",
            isTRUE(all.equal(
                unname(V_econ),
                unname(V_ref),
                tolerance = 1e-6
            )), "\n")

        cat("Maximum absolute difference:",
            max(abs(V_econ - V_ref)), "\n")
    }
}


#### 11. Ordered Model Sandwich Diagnostics ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link,
        vcov = "HC0"
    )

    S_econ = sandwich::estfun(econ)
    S_ref = sandwich::estfun(ref)

    V_econ = vcov(econ)
    V_ref = sandwich::sandwich(ref)

    cat("\nModel:", link, "\n")

    cat("Score dimensions:\n")
    print(rbind(
        EconR = dim(S_econ),
        MASS = dim(S_ref)
    ))

    cat("Score matrices equal:",
        isTRUE(all.equal(
            unname(S_econ),
            unname(S_ref),
            tolerance = 1e-6
        )), "\n")

    cat("Maximum score difference:",
        if (identical(dim(S_econ), dim(S_ref))) {
            max(abs(unname(S_econ) - unname(S_ref)))
        } else {
            NA_real_
        },
        "\n"
    )

    cat("Classical covariance matches:",
        isTRUE(all.equal(
            unname(solve(econ$Hessian)),
            unname(vcov(ref)),
            tolerance = 1e-6
        )), "\n")

    cat("HC0 standard errors:\n")

    print(data.frame(
        EconR = sqrt(diag(V_econ)),
        Sandwich = sqrt(diag(V_ref))
    ))
}


#### 12. Corrected Ordered Model HC0 Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link,
        vcov = "HC0"
    )

    V_ref = sandwich::sandwich(ref)
    V_econ = vcov(econ)

    cat("\nModel:", link, "\n")

    cat("HC0 matches sandwich:",
        isTRUE(all.equal(
            unname(V_econ),
            unname(V_ref),
            tolerance = 1e-6
        )), "\n")

    cat("Maximum absolute difference:",
        max(abs(V_econ - V_ref)), "\n")

    cat("Standard errors:\n")
    print(sqrt(diag(V_econ)))
}


#### 13. Ordered Model Bread Diagnostics ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link,
        vcov = "HC0"
    )

    n = NROW(sandwich::estfun(ref))

    B_ref = sandwich::bread(ref)
    V_classical = stats::vcov(ref)

    cat("\nModel:", link, "\n")

    cat("Bread equals n * classical covariance:",
        isTRUE(all.equal(
            unname(B_ref),
            unname(n * V_classical),
            tolerance = 1e-6
        )), "\n")

    cat("Maximum bread difference:",
        max(abs(B_ref - n * V_classical)), "\n")

    cat("Bread matrix:\n")
    print(round(B_ref, 5))

    cat("n * classical covariance:\n")
    print(round(n * V_classical, 5))

    # Reconstruct sandwich covariance from reference components
    S = sandwich::estfun(ref)

    V_reconstructed =
        B_ref %*% crossprod(S) %*% B_ref / n^2

    V_reference = sandwich::sandwich(ref)

    cat("Reconstructed sandwich matches reference:",
        isTRUE(all.equal(
            unname(V_reconstructed),
            unname(V_reference),
            tolerance = 1e-6
        )), "\n")
}


#### 14. Ordered Model Weight Normalization ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        link = link,
        vcov = "HC0"
    )

    S = sandwich::estfun(ref)

    n_rows = nrow(S)
    n_weighted = sum(econ$weights)

    B = sandwich::bread(ref)
    V_classical = stats::vcov(ref)

    cat("\nModel:", link, "\n")

    cat("Number of rows:", n_rows, "\n")
    cat("Sum of weights:", n_weighted, "\n")

    cat("Bread equals sum(weights) * vcov:",
        isTRUE(all.equal(
            unname(B),
            unname(n_weighted * V_classical),
            tolerance = 1e-6
        )), "\n")

    V_econ = vcov(econ)

    V_weighted = B %*%
        crossprod(S) %*%
        B / n_weighted^2

    cat("EconR matches weight-normalized sandwich:",
        isTRUE(all.equal(
            unname(V_econ),
            unname(V_weighted),
            tolerance = 1e-6
        )), "\n")

    cat("Maximum absolute difference:",
        max(abs(V_econ - V_weighted)), "\n")
}


#### 15. Ordered Model Unweighted HC0 Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df,
        link = link,
        vcov = "HC0"
    )

    V_ref = sandwich::sandwich(ref)
    V_econ = vcov(econ)

    cat("\nModel:", link, "\n")

    cat("HC0 matches sandwich:",
        isTRUE(all.equal(
            unname(V_econ),
            unname(V_ref),
            tolerance = 1e-6
        )), "\n")

    cat("Maximum absolute difference:",
        max(abs(V_econ - V_ref)), "\n")
}


#### 16. Ordered Model HC1 Validation ####

for (link in c("logit", "probit")) {

    for (weighted in c(FALSE, TRUE)) {

        args = list(
            formula = y ~ x1 + x2 + x3,
            data = if (weighted) df_test else df,
            link = link
        )

        if (weighted) {
            args$weights = df_test$w
            args$subset = df_test$keep
        }

        model_hc0 = do.call(
            econ_ordered,
            c(args, list(vcov = "HC0"))
        )

        model_hc1 = do.call(
            econ_ordered,
            c(args, list(vcov = "HC1"))
        )

        V0 = vcov(model_hc0)
        V1 = vcov(model_hc1)

        n = length(model_hc1$used)
        k = length(coef(model_hc1)) + length(model_hc1$zeta)

        V_expected = V0 * n / (n - k)

        cat("\nModel:", link,
            "| Weighted:", weighted, "\n")

        cat("HC1 correction:", n / (n - k), "\n")

        cat("HC1 matches expected:",
            isTRUE(all.equal(
                unname(V1),
                unname(V_expected),
                tolerance = 1e-10
            )), "\n")

        cat("Standard errors finite:",
            all(is.finite(sqrt(diag(V1)))), "\n")
    }
}


#### 17. Independent Ordered Model HC1 Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df,
        link = link,
        vcov = "HC1"
    )

    n = nrow(sandwich::estfun(ref))
    k = length(coef(ref)) + length(ref$zeta)

    V_ref = sandwich::sandwich(ref) * n / (n - k)
    V_econ = vcov(econ)

    cat("\nModel:", link, "\n")

    cat("HC1 matches independent reference:",
        isTRUE(all.equal(
            unname(V_econ),
            unname(V_ref),
            tolerance = 1e-10
        )), "\n")

    cat("Maximum absolute difference:",
        max(abs(V_econ - V_ref)), "\n")
}


#### 18. Ordered Model Cluster-Robust Reference ####

set.seed(123)

df_cluster = df
df_cluster$cluster_id = sample(
    1:50,
    nrow(df_cluster),
    replace = TRUE
)

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_cluster,
        method = method,
        Hess = TRUE
    )

    V_CR0 = sandwich::vcovCL(
        ref,
        cluster = df_cluster$cluster_id,
        type = "HC0",
        cadjust = FALSE
    )

    V_CR1 = sandwich::vcovCL(
        ref,
        cluster = df_cluster$cluster_id,
        type = "HC1",
        cadjust = FALSE
    )

    cat("\nModel:", link, "\n")
    cat("Number of clusters:",
        length(unique(df_cluster$cluster_id)), "\n")

    cat("CR0 standard errors:\n")
    print(sqrt(diag(V_CR0)))

    cat("CR1 correction verified:",
        isTRUE(all.equal(
            unname(V_CR1),
            unname(V_CR0 * nrow(df_cluster) / (nrow(df_cluster) - 5)),
            tolerance = 1e-10
        )), "\n")
}


#### 19. Ordered Model CR1 Correction Diagnostics ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_cluster,
        method = method,
        Hess = TRUE
    )

    V0 = sandwich::vcovCL(
        ref,
        cluster = df_cluster$cluster_id,
        type = "HC0",
        cadjust = FALSE
    )

    V1 = sandwich::vcovCL(
        ref,
        cluster = df_cluster$cluster_id,
        type = "HC1",
        cadjust = FALSE
    )

    cat("\nModel:", link, "\n")

    cat("nobs:", stats::nobs(ref), "\n")
    cat("Residual df:", stats::df.residual(ref), "\n")

    cat("Observed CR1/CR0 ratio:",
        V1[1, 1] / V0[1, 1], "\n")

    cat("Expected using residual df:",
        stats::nobs(ref) / stats::df.residual(ref), "\n")

    cat("Expected using 5 parameters:",
        nrow(df_cluster) / (nrow(df_cluster) - 5), "\n")
}


#### 20. Ordered Model CR1 Correction Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_cluster,
        method = method,
        Hess = TRUE
    )

    V0 = sandwich::vcovCL(
        ref,
        cluster = df_cluster$cluster_id,
        type = "HC0",
        cadjust = FALSE
    )

    V1 = sandwich::vcovCL(
        ref,
        cluster = df_cluster$cluster_id,
        type = "HC1",
        cadjust = FALSE
    )

    n = stats::nobs(ref)
    k = n - stats::df.residual(ref)

    V_expected = V0 * (n - 1) / (n - k)

    cat("\nModel:", link, "\n")

    cat("CR1 matches reference:",
        isTRUE(all.equal(
            unname(V1),
            unname(V_expected),
            tolerance = 1e-10
        )), "\n")

    cat("Maximum absolute difference:",
        max(abs(V1 - V_expected)), "\n")
}


#### 21. Investigate HC2 and HC3 Support ####

# Check available arguments
args(rms::robcov)

# Check sandwich methods for ordered models
methods("hatvalues")
methods("estfun")
methods("bread")

# Check whether HC2/HC3 are supported directly
getAnywhere("vcovHC.default")


#### 22. Ordered Model HC2/HC3 Compatibility ####

ref_polr = MASS::polr(
    y ~ x1 + x2 + x3,
    data = df,
    method = "logistic",
    Hess = TRUE
)

ref_clm = ordinal::clm(
    y ~ x1 + x2 + x3,
    data = df,
    link = "logit",
    Hess = TRUE
)

for (model_name in c("polr", "clm")) {

    model = if (model_name == "polr") ref_polr else ref_clm

    cat("\nModel:", model_name, "\n")

    for (type in c("HC0", "HC2", "HC3")) {

        result = tryCatch(
            sandwich::vcovHC(model, type = type),
            error = function(e) e
        )

        if (inherits(result, "error")) {
            cat(type, "ERROR:", conditionMessage(result), "\n")
        } else {
            cat(type, "SUCCESS\n")
            print(sqrt(diag(result)))
        }
    }
}


#### 23. Ordered Model CR0 Validation ####

devtools::load_all()

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_cluster,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_cluster,
        cluster = cluster_id,
        link = link,
        vcov = "CR0"
    )

    V_ref = sandwich::vcovCL(
        ref,
        cluster = df_cluster$cluster_id,
        type = "HC0",
        cadjust = FALSE
    )

    V_econ = vcov(econ)

    cat("\nModel:", link, "\n")

    cat("CR0 matches sandwich:",
        isTRUE(all.equal(
            unname(V_econ),
            unname(V_ref),
            tolerance = 1e-8
        )), "\n")

    cat("Maximum absolute difference:",
        max(abs(V_econ - V_ref)), "\n")

    cat("Number of clusters:",
        length(unique(econ$cluster)), "\n")
}


#### 24. Ordered Model CR1 Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_cluster,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_cluster,
        cluster = cluster_id,
        link = link,
        vcov = "CR1"
    )

    V_ref = sandwich::vcovCL(
        ref,
        cluster = df_cluster$cluster_id,
        type = "HC1",
        cadjust = FALSE
    )

    V_econ = vcov(econ)

    cat("\nModel:", link, "\n")

    cat("CR1 matches sandwich:",
        isTRUE(all.equal(
            unname(V_econ),
            unname(V_ref),
            tolerance = 1e-8
        )), "\n")

    cat("Maximum absolute difference:",
        max(abs(V_econ - V_ref)), "\n")
}


#### 25. Weighted and Subsetted Cluster Validation ####

set.seed(456)

df_test$cluster_id = sample(
    1:50,
    nrow(df_test),
    replace = TRUE
)

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        method = method,
        Hess = TRUE
    )

    for (type in c("CR0", "CR1")) {

        econ = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df_test,
            weights = w,
            subset = keep,
            cluster = cluster_id,
            link = link,
            vcov = type
        )

        V_ref = sandwich::vcovCL(
            ref,
            cluster = econ$cluster,
            type = if (type == "CR0") "HC0" else "HC1",
            cadjust = FALSE
        )

        V_econ = vcov(econ)

        cat("\nModel:", link, "| Covariance:", type, "\n")

        cat("Matches sandwich:",
            isTRUE(all.equal(
                unname(V_econ),
                unname(V_ref),
                tolerance = 1e-8
            )), "\n")

        cat("Maximum difference:",
            max(abs(V_econ - V_ref)), "\n")

        cat("Estimation rows:", length(econ$used), "\n")
        cat("Number of clusters:", length(unique(econ$cluster)), "\n")
    }
}


#### 26. Weighted Cluster Covariance Diagnostics ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        cluster = cluster_id,
        link = link,
        vcov = "CR0"
    )

    S = sandwich::estfun(ref)

    n = nrow(S)
    k = ncol(S)
    n_weighted = sum(econ$weights)

    cluster_scores = rowsum(
        S,
        group = econ$cluster,
        reorder = TRUE
    )

    V_classical = stats::vcov(ref)

    V_manual = V_classical %*%
        crossprod(cluster_scores) %*%
        V_classical

    V_econ = vcov(econ)

    V_ref = sandwich::vcovCL(
        ref,
        cluster = econ$cluster,
        type = "HC0",
        cadjust = FALSE
    )

    cat("\nModel:", link, "\n")

    cat("Cluster alignment:",
        identical(
            econ$cluster,
            df_test$cluster_id[econ$used]
        ), "\n")

    cat("EconR matches manual CR0:",
        isTRUE(all.equal(
            unname(V_econ),
            unname(V_manual),
            tolerance = 1e-8
        )), "\n")

    cat("Reference / EconR covariance ratio:",
        V_ref[1, 1] / V_econ[1, 1], "\n")

    cat("Expected normalization ratio:",
        (n_weighted / n)^2, "\n")

    cat("Reference matches scaled EconR:",
        isTRUE(all.equal(
            unname(V_ref),
            unname(V_econ * (n_weighted / n)^2),
            tolerance = 1e-8
        )), "\n")
}


#### 27. Weighted CR1 Validation ####

for (link in c("logit", "probit")) {

    model_cr0 = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        cluster = cluster_id,
        link = link,
        vcov = "CR0"
    )

    model_cr1 = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_test,
        weights = w,
        subset = keep,
        cluster = cluster_id,
        link = link,
        vcov = "CR1"
    )

    n = length(model_cr1$used)
    k = length(coef(model_cr1)) + length(model_cr1$zeta)

    V_expected = vcov(model_cr0) * (n - 1) / (n - k)

    cat("\nModel:", link, "\n")

    cat("CR1 correction:", (n - 1) / (n - k), "\n")

    cat("CR1 matches expected:",
        isTRUE(all.equal(
            unname(vcov(model_cr1)),
            unname(V_expected),
            tolerance = 1e-10
        )), "\n")

    cat("Maximum absolute difference:",
        max(abs(vcov(model_cr1) - V_expected)), "\n")
}


#### 28. Ordered Model Summary Diagnostics ####

for (link in c("logit", "probit")) {

    for (type in c("classical", "HC0", "HC1", "CR0", "CR1")) {

        model = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df_cluster,
            link = link,
            cluster = cluster_id,
            vcov = type
        )

        s = summary(model)

        # Standard errors reported by summary()
        se_summary = s$coefficients[, "Std. Error"]

        # Standard errors from EconR's covariance estimator
        se_vcov = sqrt(diag(vcov(model)))

        cat("\nModel:", link, "| Covariance:", type, "\n")

        cat("Summary matches EconR covariance:",
            isTRUE(all.equal(
                unname(se_summary),
                unname(se_vcov),
                tolerance = 1e-8
            )), "\n")

        cat("Maximum SE difference:",
            max(abs(se_summary - se_vcov)), "\n")
    }
}


#### 29. Ordered Model Prediction Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df,
        method = method,
        Hess = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df,
        link = link
    )

    # Predicted category probabilities
    prob_ref = predict(ref, type = "probs")
    prob_econ = predict(econ, type = "probs")

    # Predicted outcome categories
    class_ref = predict(ref, type = "class")
    class_econ = predict(econ, type = "class")

    cat("\nModel:", link, "\n")

    cat("Predicted probabilities match:",
        isTRUE(all.equal(
            unname(prob_econ),
            unname(prob_ref),
            tolerance = 1e-8
        )), "\n")

    cat("Predicted categories match:",
        identical(class_econ, class_ref), "\n")

    cat("Probability rows sum to one:",
        all(abs(rowSums(prob_econ) - 1) < 1e-10), "\n")

    cat("Prediction dimensions:",
        paste(dim(prob_econ), collapse = " x "), "\n")
}


#### 30. Ordered Model Methods Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df,
        method = method,
        Hess = TRUE,
        model = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df,
        link = link
    )

    cat("\nModel:", link, "\n")

    # Test each method without interrupting the loop
    checks = list(
        logLik = function(m) stats::logLik(m),
        nobs = function(m) stats::nobs(m),
        fitted = function(m) stats::fitted(m),
        residuals = function(m) stats::residuals(m),
        weights = function(m) stats::weights(m),
        model.matrix = function(m) stats::model.matrix(m)
    )

    for (method_name in names(checks)) {

        result = tryCatch({

            value_ref = checks[[method_name]](ref)
            value_econ = checks[[method_name]](econ)

            isTRUE(all.equal(
                value_econ,
                value_ref,
                tolerance = 1e-8
            ))

        }, error = function(e) {
            paste("ERROR:", conditionMessage(e))
        })

        cat(method_name, ":", result, "\n")
    }
}


#### 31. Ordered Model Weights Diagnostics ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df,
        method = method,
        Hess = TRUE,
        model = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df,
        link = link
    )

    cat("\nModel:", link, "\n")

    cat("Reference weights:\n")
    print(head(stats::weights(ref)))

    cat("EconR weights:\n")
    print(head(stats::weights(econ)))

    cat("Reference weights NULL:",
        is.null(stats::weights(ref)), "\n")

    cat("EconR weights NULL:",
        is.null(stats::weights(econ)), "\n")

    cat("Stored EconR weights match expected:",
        isTRUE(all.equal(
            unname(econ$weights),
            rep(1, nrow(df))
        )), "\n")
}


#### 32. Ordered Model Tidy Validation ####

for (link in c("logit", "probit")) {

    for (type in c("classical", "HC0", "HC1", "CR0", "CR1")) {

        model = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df_cluster,
            link = link,
            cluster = cluster_id,
            vcov = type
        )

        result = tryCatch(
            broom::tidy(model),
            error = function(e) e
        )

        cat("\nModel:", link, "| Covariance:", type, "\n")

        if (inherits(result, "error")) {

            cat("ERROR:", conditionMessage(result), "\n")

        } else {

            se_tidy = result$std.error
            se_vcov = sqrt(diag(vcov(model)))

            cat("Tidy standard errors match EconR:",
                isTRUE(all.equal(
                    unname(se_tidy),
                    unname(se_vcov),
                    tolerance = 1e-8
                )), "\n")

            cat("Number of coefficients:", nrow(result), "\n")
        }
    }
}


#### 33. Ordered Model Glance and Augment Validation ####

for (link in c("logit", "probit")) {

    method = if (link == "logit") "logistic" else "probit"

    ref = MASS::polr(
        y ~ x1 + x2 + x3,
        data = df,
        method = method,
        Hess = TRUE,
        model = TRUE
    )

    econ = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df,
        link = link,
        vcov = "HC1"
    )

    cat("\nModel:", link, "\n")

    # Model-level statistics
    glance_ref = tryCatch(
        broom::glance(ref),
        error = function(e) e
    )

    glance_econ = tryCatch(
        broom::glance(econ),
        error = function(e) e
    )

    if (inherits(glance_econ, "error")) {

        cat("Glance ERROR:",
            conditionMessage(glance_econ), "\n")

    } else if (inherits(glance_ref, "error")) {

        cat("Reference glance ERROR:",
            conditionMessage(glance_ref), "\n")

    } else {

        cat("Glance matches MASS:",
            isTRUE(all.equal(
                glance_econ,
                glance_ref,
                tolerance = 1e-8
            )), "\n")

        print(glance_econ)
    }

    # Observation-level output
    augment_econ = tryCatch(
        broom::augment(econ),
        error = function(e) e
    )

    if (inherits(augment_econ, "error")) {

        cat("Augment ERROR:",
            conditionMessage(augment_econ), "\n")

    } else {

        cat("Augment rows:", nrow(augment_econ), "\n")
        cat("Augment columns:\n")
        print(names(augment_econ))
    }
}


#### 34. Ordered Model Modelsummary Validation ####

models = list(
    "Ordered Logit" = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_cluster,
        link = "logit",
        cluster = cluster_id,
        vcov = "CR1"
    ),

    "Ordered Probit" = econ_ordered(
        y ~ x1 + x2 + x3,
        data = df_cluster,
        link = "probit",
        cluster = cluster_id,
        vcov = "CR1"
    )
)

# Display regression table
modelsummary::modelsummary(
    models,
    estimate = "{estimate}",
    statistic = "({std.error})",
    output = "markdown"
)

# Verify the standard errors independently
for (name in names(models)) {

    model = models[[name]]
    tidy = broom::tidy(model)

    cat("\nModel:", name, "\n")

    cat("Modelsummary input SE matches EconR:",
        isTRUE(all.equal(
            unname(tidy$std.error),
            unname(sqrt(diag(vcov(model)))),
            tolerance = 1e-8
        )), "\n")
}


#### 35. Test clubSandwich support for ordered models ####

set.seed(123)

n = 500

df = data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    cluster = rep(seq_len(25), each = 20)
)

latent = 0.8 * df$x1 - 0.4 * df$x2 + rlogis(n)

df$y = ordered(
    cut(
        latent,
        breaks = c(-Inf, -1, 0, 1, Inf),
        labels = c("Low", "Medium", "High", "Very high")
    )
)

reference = MASS::polr(
    y ~ x1 + x2,
    data = df,
    method = "logistic",
    Hess = TRUE
)

for (type in c("CR0", "CR1", "CR1p", "CR1S", "CR2", "CR3")) {

    cat("\n", type, ": ", sep = "")

    tryCatch(
        {
            V = clubSandwich::vcovCR(
                reference,
                cluster = df$cluster,
                type = type
            )

            cat("SUPPORTED\n")
            print(dim(V))
        },
        error = function(e) {
            cat("NOT SUPPORTED: ", conditionMessage(e), "\n")
        }
    )
}


#### 36. Validate ordered-model scores ####

library(numDeriv)

model = econ_ordered(
    y ~ x1 + x2,
    data = df,
    link = "logit",
    vcov = "classical"
)

S = sandwich::estfun(model)

X = stats::model.matrix(model)
X = X[, colnames(X) != "(Intercept)", drop = FALSE]

y = as.integer(stats::model.response(model$model))

K = length(stats::coef(model))
J = length(model$zeta) + 1L

theta = c(stats::coef(model), model$zeta)

loglik_i = function(par, i) {

    beta = par[seq_len(K)]
    zeta = par[K + seq_len(J - 1L)]

    eta = drop(X[i, , drop = FALSE] %*% beta)

    upper = if (y[i] == J) Inf else zeta[y[i]]
    lower = if (y[i] == 1L) -Inf else zeta[y[i] - 1L]

    p = plogis(upper - eta) - plogis(lower - eta)

    log(p)
}

indices = c(1, 2, 3, 10, 25, 50)

comparison = do.call(
    rbind,
    lapply(indices, function(i) {

        numerical = numDeriv::grad(
            func = function(par) loglik_i(par, i),
            x = theta
        )

        data.frame(
            observation = i,
            parameter = names(theta),
            analytical = as.numeric(S[i, ]),
            numerical = as.numeric(numerical),
            difference = as.numeric(S[i, ]) - as.numeric(numerical)
        )
    })
)

print(comparison)

cat(
    "\nMaximum absolute difference:",
    max(abs(comparison$difference)),
    "\n"
)






















#### Dev Tool Testing ####
devtools::load_all()
devtools::test()

devtools::document()

devtools::check()
