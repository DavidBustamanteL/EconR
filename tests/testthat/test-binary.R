#### 1. Binary Model Estimation ####

test_that("Logit estimation matches glm", {

    econ = econ_binary(
        am ~ wt + hp,
        data = mtcars,
        link = "logit"
    )

    reference = glm(
        am ~ wt + hp,
        data = mtcars,
        family = binomial(link = "logit")
    )

    expect_equal(
        coef(econ),
        coef(reference),
        tolerance = 1e-10
    )

    expect_equal(
        vcov(econ),
        vcov(reference),
        tolerance = 1e-10
    )

    expect_equal(
        econ$fitted.values,
        fitted(reference),
        tolerance = 1e-10
    )
})

#### 2. Clustered Covariance ####

test_that("Binary clustered covariance matches reference estimators", {

    glm_ref = glm(
        am ~ wt + hp,
        data = mtcars,
        family = binomial(link = "logit")
    )

    for (type in c("CR0", "CR2", "stata")) {

        econ = econ_binary(
            am ~ wt + hp,
            data = mtcars,
            link = "logit",
            vcov = type,
            cluster = cyl
        )

        V_ref = switch(
            type,
            CR0 = clubSandwich::vcovCR(
                glm_ref,
                cluster = mtcars$cyl,
                type = "CR0"
            ),
            CR2 = clubSandwich::vcovCR(
                glm_ref,
                cluster = mtcars$cyl,
                type = "CR2"
            ),
            stata = sandwich::vcovCL(
                glm_ref,
                cluster = mtcars$cyl,
                type = "HC1"
            )
        )

        tolerance = if (type == "CR0") 1e-6 else 1e-10

        expect_equal(
            unclass(vcov(econ)),
            unclass(as.matrix(V_ref)),
            tolerance = tolerance
        )

        expect_equal(
            unname(summary(econ)$coefficients[, 2]),
            unname(sqrt(diag(vcov(econ)))),
            tolerance = 1e-10
        )
    }
})


#### 3. Probit Model Estimation ####

test_that("Probit estimation matches glm", {

    expect_warning(
        econ <- econ_binary(
            am ~ wt + hp,
            data = mtcars,
            link = "probit"
        ),
        "fitted probabilities numerically 0 or 1 occurred"
    )

    expect_warning(
        reference <- glm(
            am ~ wt + hp,
            data = mtcars,
            family = binomial(link = "probit")
        ),
        "fitted probabilities numerically 0 or 1 occurred"
    )

    expect_equal(
        coef(econ),
        coef(reference),
        tolerance = 1e-10
    )

    expect_equal(
        vcov(econ),
        vcov(reference),
        tolerance = 1e-10
    )

    expect_equal(
        econ$fitted.values,
        fitted(reference),
        tolerance = 1e-10
    )
})

#### 4. Binary Model Predictions ####

test_that("Binary predictions match glm", {

    for (link in c("logit", "probit")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link
            )
        )

        reference = suppressWarnings(
            glm(
                am ~ wt + hp,
                data = mtcars,
                family = binomial(link = link)
            )
        )

        #### 4.1. In-sample predictions ####

        expect_equal(
            predict(econ, type = "link"),
            predict(reference, type = "link"),
            tolerance = 1e-10
        )

        expect_equal(
            predict(econ, type = "response"),
            predict(reference, type = "response"),
            tolerance = 1e-10
        )

        #### 4.2. Out-of-sample predictions ####

        newdata = data.frame(
            wt = c(2.5, 3.0, 3.5),
            hp = c(100, 150, 200)
        )

        expect_equal(
            predict(econ, newdata = newdata, type = "link"),
            predict(reference, newdata = newdata, type = "link"),
            tolerance = 1e-10
        )

        expect_equal(
            predict(econ, newdata = newdata, type = "response"),
            predict(reference, newdata = newdata, type = "response"),
            tolerance = 1e-10
        )
    }
})

#### 5. Log-Likelihood and Information Criteria ####

test_that("Binary likelihood and information criteria match glm", {

    for (link in c("logit", "probit")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link
            )
        )

        reference = suppressWarnings(
            glm(
                am ~ wt + hp,
                data = mtcars,
                family = binomial(link = link)
            )
        )

        #### 5.1. Log-likelihood ####

        expect_equal(
            as.numeric(logLik(econ)),
            as.numeric(logLik(reference)),
            tolerance = 1e-10
        )

        #### 5.2. AIC ####

        expect_equal(
            AIC(econ),
            AIC(reference),
            tolerance = 1e-10
        )

        #### 5.3. BIC ####

        expect_equal(
            BIC(econ),
            BIC(reference),
            tolerance = 1e-10
        )
    }
})

#### 6. Weighted Binary Estimation ####

test_that("Weighted binary estimation matches glm", {

    df = mtcars
    df$w = df$disp / mean(df$disp)

    for (link in c("logit", "probit")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = df,
                weights = w,
                link = link
            )
        )

        reference = suppressWarnings(
            glm(
                am ~ wt + hp,
                data = df,
                weights = w,
                family = binomial(link = link)
            )
        )

        #### 6.1. Coefficients ####

        expect_equal(
            coef(econ),
            coef(reference),
            tolerance = 1e-10
        )

        #### 6.2. Classical covariance ####

        expect_equal(
            vcov(econ),
            vcov(reference),
            tolerance = 1e-10
        )

        #### 6.3. Fitted probabilities ####

        expect_equal(
            econ$fitted.values,
            fitted(reference),
            tolerance = 1e-10
        )

        #### 6.4. Observation weights ####

        expect_equal(
            as.numeric(econ$weights),
            as.numeric(model.weights(model.frame(reference))),
            tolerance = 1e-10
        )
    }
})

#### 7. Binary Sample Handling ####

test_that("Binary sample selection matches glm", {

    for (link in c("logit", "probit")) {

        df = mtcars

        #### 7.1. Subset ####

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = df,
                subset = mpg > 18,
                link = link
            )
        )

        reference = suppressWarnings(
            glm(
                am ~ wt + hp,
                data = df,
                subset = mpg > 18,
                family = binomial(link = link)
            )
        )

        expect_equal(coef(econ), coef(reference), tolerance = 1e-10)
        expect_equal(vcov(econ), vcov(reference), tolerance = 1e-10)
        expect_equal(nobs(econ), nobs(reference))

        #### 7.2. Missing covariates ####

        df_na = df
        df_na$hp[c(2, 5, 7, 12)] = NA

        econ_na = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = df_na,
                link = link
            )
        )

        reference_na = suppressWarnings(
            glm(
                am ~ wt + hp,
                data = df_na,
                family = binomial(link = link)
            )
        )

        expect_equal(coef(econ_na), coef(reference_na), tolerance = 1e-10)
        expect_equal(vcov(econ_na), vcov(reference_na), tolerance = 1e-10)
        expect_equal(nobs(econ_na), nobs(reference_na))

        #### 7.3. Missing weights ####

        df_w = df
        df_w$w = rep(1, nrow(df_w))
        df_w$w[c(3, 9)] = NA

        econ_w = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = df_w,
                weights = w,
                link = link
            )
        )

        reference_w = suppressWarnings(
            glm(
                am ~ wt + hp,
                data = df_w,
                weights = w,
                family = binomial(link = link)
            )
        )

        expect_equal(coef(econ_w), coef(reference_w), tolerance = 1e-10)
        expect_equal(vcov(econ_w), vcov(reference_w), tolerance = 1e-10)
        expect_equal(nobs(econ_w), nobs(reference_w))
    }
})

#### 8. Heteroskedasticity-Robust Covariance ####

test_that("Binary HC0-HC3 covariance matches sandwich", {

    for (link in c("logit", "probit")) {

        reference = suppressWarnings(
            glm(
                am ~ wt + hp,
                data = mtcars,
                family = binomial(link = link)
            )
        )

        for (type in c("HC0", "HC1", "HC2", "HC3")) {

            econ = suppressWarnings(
                econ_binary(
                    am ~ wt + hp,
                    data = mtcars,
                    link = link,
                    vcov = type
                )
            )

            #### 8.1. Covariance matrix ####

            expect_equal(
                vcov(econ),
                sandwich::vcovHC(reference, type = type),
                tolerance = 1e-10
            )

            #### 8.2. Summary standard errors ####

            expect_equal(
                unname(summary(econ)$coefficients[, 2]),
                unname(sqrt(diag(vcov(econ)))),
                tolerance = 1e-10
            )
        }
    }
})

#### 9. Sandwich Compatibility Methods ####

test_that("Binary sandwich methods match glm", {

    for (link in c("logit", "probit")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link
            )
        )

        reference = suppressWarnings(
            glm(
                am ~ wt + hp,
                data = mtcars,
                family = binomial(link = link)
            )
        )

        #### 9.1. Bread matrix ####

        expect_equal(
            sandwich::bread(econ),
            sandwich::bread(reference),
            tolerance = 1e-10
        )

        #### 9.2. Estimating functions ####

        expect_equal(
            sandwich::estfun(econ),
            sandwich::estfun(reference),
            tolerance = 1e-10
        )

        #### 9.3. Leverage values ####

        expect_equal(
            hatvalues(econ),
            hatvalues(reference),
            tolerance = 1e-10
        )

        #### 9.4. Model matrix ####

        expect_equal(
            model.matrix(econ),
            model.matrix(reference),
            tolerance = 1e-10
        )
    }
})

#### 10. Clustered Covariance Sample Alignment ####

test_that("Binary CR2 respects subsets and missing observations", {

    df = mtcars
    df$hp[c(2, 7)] = NA

    econ = suppressWarnings(
        econ_binary(
            am ~ wt + hp,
            data = df,
            subset = mpg > 18,
            link = "logit",
            vcov = "CR2",
            cluster = cyl
        )
    )

    reference = suppressWarnings(
        glm(
            am ~ wt + hp,
            data = econ$model,
            family = binomial(link = "logit")
        )
    )

    V_reference = clubSandwich::vcovCR(
        reference,
        cluster = econ$cluster,
        type = "CR2"
    )

    #### 10.1. Sample alignment ####

    expect_equal(
        length(econ$cluster),
        nobs(econ)
    )

    #### 10.2. Coefficients ####

    expect_equal(
        coef(econ),
        coef(reference),
        tolerance = 1e-10
    )

    #### 10.3. CR2 covariance ####

    expect_equal(
        unclass(vcov(econ)),
        unclass(as.matrix(V_reference)),
        tolerance = 1e-10
    )
})

#### 11. Probit Covariance Estimation ####

test_that("Probit covariance matches reference estimators", {

    reference = suppressWarnings(
        glm(
            am ~ wt + hp,
            data = mtcars,
            family = binomial(link = "probit")
        )
    )

    #### 11.1. Classical covariance ####

    econ = suppressWarnings(
        econ_binary(
            am ~ wt + hp,
            data = mtcars,
            link = "probit"
        )
    )

    expect_equal(
        vcov(econ),
        vcov(reference),
        tolerance = 1e-10
    )

    #### 11.2. Clustered covariance ####

    for (type in c("CR0", "CR2", "stata")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = "probit",
                vcov = type,
                cluster = cyl
            )
        )

        V_reference = switch(
            type,
            CR0 = sandwich::vcovCL(
              reference,
              cluster = mtcars$cyl,
              type = "HC",
              cadjust = FALSE
            ),
            CR2 = clubSandwich::vcovCR(
                reference,
                cluster = mtcars$cyl,
                type = "CR2"
            ),
            stata = sandwich::vcovCL(
                reference,
                cluster = mtcars$cyl,
                type = "HC1"
            )
        )

        expect_equal(
            unclass(suppressWarnings(vcov(econ))),
            unclass(as.matrix(V_reference)),
            tolerance = 1e-10
        )
    }
})

#### 12. Probit Missing-Data Handling ####

test_that("Probit handles missing data and cluster alignment", {

    df = mtcars
    df$hp[c(2, 7)] = NA
    df$w = rep(1, nrow(df))
    df$w[c(3, 9)] = NA

    econ = suppressWarnings(
        econ_binary(
            am ~ wt + hp,
            data = df,
            weights = w,
            subset = mpg > 18,
            link = "probit",
            vcov = "CR2",
            cluster = cyl
        )
    )

    reference = suppressWarnings(
        glm(
            am ~ wt + hp,
            data = df,
            weights = w,
            subset = mpg > 18,
            family = binomial(link = "probit")
        )
    )

    #### 12.1. Estimation ####

    expect_equal(coef(econ), coef(reference), tolerance = 1e-10)
    expect_equal(nobs(econ), nobs(reference))

    #### 12.2. Sample alignment ####

    expect_equal(length(econ$cluster), nobs(econ))
    expect_equal(length(econ$weights), nobs(econ))

    #### 12.3. Clustered covariance ####

    V_reference = suppressWarnings(
        clubSandwich::vcovCR(
            reference,
            cluster = econ$cluster,
            type = "CR2"
        )
    )

    expect_equal(
        unclass(suppressWarnings(vcov(econ))),
        unclass(as.matrix(V_reference)),
        tolerance = 1e-10
    )
})

#### 13. Binary Model Methods ####

test_that("Binary model methods match glm", {

    for (link in c("logit", "probit")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link
            )
        )

        reference = suppressWarnings(
            glm(
                am ~ wt + hp,
                data = mtcars,
                family = binomial(link = link)
            )
        )

        #### 13.1. Coefficients and covariance ####

        expect_equal(coef(econ), coef(reference), tolerance = 1e-10)
        expect_equal(vcov(econ), vcov(reference), tolerance = 1e-10)

        #### 13.2. Model dimensions ####

        expect_equal(nobs(econ), nobs(reference))

        #### 13.3. Likelihood methods ####

        expect_equal(
            attr(logLik(econ), "df"),
            attr(logLik(reference), "df")
        )

        expect_equal(
            attr(logLik(econ), "nobs"),
            attr(logLik(reference), "nobs")
        )

        #### 13.4. Prediction interface ####

        expect_equal(
            predict(econ, type = "response"),
            predict(reference, type = "response"),
            tolerance = 1e-10
        )
    }
})

#### 14. Summary and Inference ####

test_that("Binary summary and inference are correct", {

    for (link in c("logit", "probit")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link,
                vcov = "HC3"
            )
        )

        S = summary(econ)$coefficients

        #### 14.1. Coefficient estimates ####

        expect_equal(
            unname(S[, 1]),
            unname(coef(econ)),
            tolerance = 1e-10
        )

        #### 14.2. Standard errors ####

        SE = sqrt(diag(vcov(econ)))

        expect_equal(
            unname(S[, 2]),
            unname(SE),
            tolerance = 1e-10
        )

        #### 14.3. Z-statistics ####

        z = coef(econ) / SE

        expect_equal(
            unname(S[, 3]),
            unname(z),
            tolerance = 1e-10
        )

        #### 14.4. P-values ####

        p = 2 * stats::pnorm(abs(z), lower.tail = FALSE)

        expect_equal(
            unname(S[, 4]),
            unname(p),
            tolerance = 1e-10
        )

        #### 14.5. Printed summary ####

        expect_output(
            print(summary(econ)),
            "Coefficients"
        )
    }
})

#### 15. Broom and Modelsummary Compatibility ####

test_that("Binary models support broom and modelsummary", {

    skip_if_not_installed("broom")
    skip_if_not_installed("modelsummary")

    for (link in c("logit", "probit")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link,
                vcov = "HC3"
            )
        )

        #### 15.1. Tidy output ####

        T = broom::tidy(econ)

        expect_s3_class(T, "data.frame")
        expect_equal(T$estimate, unname(coef(econ)))
        expect_equal(
            T$std.error,
            unname(sqrt(diag(vcov(econ)))),
            tolerance = 1e-10
        )

        #### 15.2. Glance output ####

        G = broom::glance(econ)

        expect_equal(G$nobs, nobs(econ))
        expect_equal(G$aic, AIC(econ))
        expect_equal(G$bic, BIC(econ))
        expect_equal(G$logLik, as.numeric(logLik(econ)))

        #### 15.3. Modelsummary output ####

        M = modelsummary::modelsummary(
            list("EconR" = econ),
            output = "data.frame"
        )

        expect_s3_class(M, "data.frame")
        expect_true(all(
            c("(Intercept)", "wt", "hp") %in% M$term
        ))
    }
})


#### 16. Interaction and Quadratic Terms ####

set.seed(123)

n = 600

wages = data.frame(
    educ = sample(8:20, n, replace = TRUE),
    exper = sample(0:40, n, replace = TRUE),
    nonwhite = rbinom(n, 1, 0.25),
    female = rbinom(n, 1, 0.50),
    wage = runif(n, 5, 35)
)

eta = with(
    wages,
    -2 + 0.08 * educ - 0.3 * female +
        0.015 * educ * female +
        0.12 * exper - 0.002 * exper^2 +
        0.03 * wage - 0.4 * nonwhite
)

wages$married = rbinom(n, 1, plogis(eta))

    formulas = list(
        interaction = married ~ educ * female + exper + nonwhite + wage,
        quadratic = married ~ educ + exper + I(exper^2) + nonwhite + female + wage,
        combined = married ~ educ * female + exper + I(exper^2) + nonwhite + wage
    )

    for (f in formulas) {

        ref = stats::glm(
            f,
            data = wages,
            family = stats::binomial(link = "logit")
        )

        econ = econ_binary(
            f,
            data = wages,
            link = "logit",
            vcov = "HC1"
        )

        expect_equal(
            coef(econ),
            coef(ref),
            tolerance = 1e-10
        )

        expect_equal(
            vcov(econ),
            sandwich::vcovHC(ref, type = "HC1"),
            tolerance = 1e-10
        )

        expect_equal(
            unname(predict(econ, type = "response")),
            unname(predict(ref, type = "response")),
            tolerance = 1e-10
        )
    }



#### 17. Binary Confidence Intervals ####

test_that("Binary confidence intervals use the selected covariance", {

    for (link in c("logit", "probit")) {

        fit = suppressWarnings(
            econ_binary(
            am ~ wt + hp,
            data = mtcars,
            link = link,
            vcov = "HC3"
        )
    )

        b = coef(fit)
        se = sqrt(diag(vcov(fit)))
        z = stats::qnorm(0.975)

        expected = cbind(
            b - z * se,
            b + z * se
        )

        colnames(expected) = c("2.5 %", "97.5 %")

        expect_equal(
            confint(fit),
            expected,
            tolerance = 1e-10
        )

        expect_equal(
            confint(fit, parm = "wt"),
            expected["wt", , drop = FALSE],
            tolerance = 1e-10
        )

        expect_error(confint(fit, level = 1))
        expect_error(confint(fit, parm = "unknown"))
    }
})


#### 18. Binary Fitted Values ####

test_that("Binary fitted values match glm", {

    for (link in c("logit", "probit")) {

        ref = suppressWarnings(
            stats::glm(
                am ~ wt + hp,
                data = mtcars,
                family = stats::binomial(link = link)
            )
        )

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link,
                vcov = "HC3"
            )
        )

        expect_equal(
            fitted(econ),
            fitted(ref),
            tolerance = 1e-10
        )

        expect_equal(
            fitted(econ),
            predict(econ, type = "response"),
            tolerance = 1e-10
        )
    }
})



#### 19. Binary Residuals ####

test_that("Binary residuals match glm", {

    for (link in c("logit", "probit")) {

        ref = suppressWarnings(
            stats::glm(
                am ~ wt + hp,
                data = mtcars,
                family = stats::binomial(link = link)
            )
        )

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link,
                vcov = "HC3"
            )
        )

        for (type in c(
            "response",
            "pearson",
            "deviance",
            "working",
            "partial"
        )) {

            expect_equal(
                unname(suppressWarnings(residuals(econ, type = type))),
                unname(suppressWarnings(residuals(ref, type = type))),
                tolerance = 1e-10
            )
        }
    }
})



#### 20. Binary Model Weights ####

test_that("Binary prior and working weights match glm", {

    data = mtcars
    data$w = data$disp / mean(data$disp)

    for (link in c("logit", "probit")) {

        ref = suppressWarnings(
            stats::glm(
                am ~ wt + hp,
                data = data,
                weights = w,
                family = stats::binomial(link = link)
            )
        )

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = data,
                weights = w,
                link = link,
                vcov = "HC1"
            )
        )

        expect_equal(
            unname(weights(econ, type = "prior")),
            unname(weights(ref, type = "prior")),
            tolerance = 1e-10
        )

        expect_equal(
            unname(weights(econ, type = "working")),
            unname(weights(ref, type = "working")),
            tolerance = 1e-10
        )
    }
})



#### 21. Binary Model Formula ####

test_that("Binary formula extraction matches glm", {

    f = am ~ wt * hp + I(wt^2)

    for (link in c("logit", "probit")) {

        ref = suppressWarnings(
            stats::glm(
                f,
                data = mtcars,
                family = stats::binomial(link = link)
            )
        )

        econ = suppressWarnings(
            econ_binary(
                f,
                data = mtcars,
                link = link,
                vcov = "HC1"
            )
        )

        expect_equal(
            formula(econ),
            formula(ref)
        )
    }
})



#### 22. Binary Analysis of Deviance ####

test_that("Single-model binary ANOVA matches glm", {

    for (link in c("logit", "probit")) {

        ref = suppressWarnings(
            stats::glm(
                am ~ wt + hp,
                data = mtcars,
                family = stats::binomial(link = link)
            )
        )

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link,
                vcov = "HC1"
            )
        )

        expect_equal(
            as.data.frame(
                suppressWarnings(anova(econ, test = "Chisq"))
            ),
            as.data.frame(
                suppressWarnings(anova(ref, test = "Chisq"))
            ),
            tolerance = 1e-10
        )
    }
})


#### 23. Nested Binary Model Comparisons ####

test_that("Nested binary ANOVA matches glm", {

    for (link in c("logit", "probit")) {

        ref1 = suppressWarnings(
            stats::glm(
                am ~ wt,
                data = mtcars,
                family = stats::binomial(link = link)
            )
        )

        ref2 = suppressWarnings(
            stats::glm(
                am ~ wt + hp,
                data = mtcars,
                family = stats::binomial(link = link)
            )
        )

        econ1 = suppressWarnings(
            econ_binary(
                am ~ wt,
                data = mtcars,
                link = link,
                vcov = "HC1"
            )
        )

        econ2 = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link,
                vcov = "HC1"
            )
        )

        expect_equal(
            as.data.frame(
                suppressWarnings(anova(econ1, econ2, test = "Chisq"))
            ),
            as.data.frame(
                suppressWarnings(anova(ref1, ref2, test = "Chisq"))
            ),
            tolerance = 1e-10
        )
    }
})


test_that("Binary ANOVA rejects incompatible estimation samples", {

    data = mtcars
    data$hp[1] = NA

    econ1 = econ_binary(
        am ~ wt,
        data = data,
        link = "logit"
    )

    econ2 = econ_binary(
        am ~ wt + hp,
        data = data,
        link = "logit"
    )

    expect_error(
        anova(econ1, econ2, test = "Chisq"),
        "identical estimation samples"
    )
})


test_that("Binary ANOVA rejects different dependent-variable values", {

    data1 = mtcars
    data2 = mtcars

    # Same sample and regressors, but a different outcome
    data2$am[1] = 1 - data2$am[1]

    econ1 = econ_binary(
        am ~ wt,
        data = data1,
        link = "logit"
    )

    econ2 = econ_binary(
        am ~ wt + hp,
        data = data2,
        link = "logit"
    )

    expect_error(
        anova(econ1, econ2, test = "Chisq"),
        "identical dependent-variable values"
    )
})



#### 24. Binary Model Augmentation ####

test_that("Binary augmentation matches fitted values and residuals", {

    for (link in c("logit", "probit")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link,
                vcov = "HC1"
            )
        )

        aug = suppressWarnings(
            broom::augment(econ, se_fit = TRUE)
        )

        expect_equal(
            aug$.fitted,
            unname(predict(econ, type = "link")),
            tolerance = 1e-10
        )

        expect_equal(
            aug$.resid,
            unname(residuals(econ, type = "deviance")),
            tolerance = 1e-10
        )

        expect_equal(
            nrow(aug),
            nobs(econ)
        )
    }
})


test_that("Binary augmentation uses EconR covariance for prediction SEs", {

    for (link in c("logit", "probit")) {

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = mtcars,
                link = link,
                vcov = "HC1"
            )
        )

        aug = suppressWarnings(
            broom::augment(econ, se_fit = TRUE)
        )

        X = model.matrix(econ)
        V = vcov(econ)

        se_manual = sqrt(pmax(
            rowSums((X %*% V) * X),
            0
        ))

        expect_equal(
            aug$.se.fit,
            unname(se_manual),
            tolerance = 1e-10
        )
    }
})


test_that("Unweighted binary augmentation has standard broom columns", {

    econ = econ_binary(
        am ~ wt + hp,
        data = mtcars,
        link = "logit",
        vcov = "HC1"
    )

    aug = broom::augment(econ, se_fit = TRUE)

    expect_named(
        aug,
        c(
            ".rownames",
            "am",
            "wt",
            "hp",
            ".fitted",
            ".se.fit",
            ".resid",
            ".hat",
            ".sigma",
            ".cooksd",
            ".std.resid"
        )
    )

    expect_false("(weights)" %in% names(aug))
})



#### 25. Weighted Binary Augmentation ####

test_that("Weighted binary augmentation matches glm", {

    df = mtcars
    df$w = df$disp / mean(df$disp)

    for (link in c("logit", "probit")) {

        ref = suppressWarnings(
            stats::glm(
                am ~ wt + hp,
                data = df,
                weights = w,
                family = stats::binomial(link = link)
            )
        )

        econ = suppressWarnings(
            econ_binary(
                am ~ wt + hp,
                data = df,
                weights = w,
                link = link,
                vcov = "HC1"
            )
        )

        aug = suppressWarnings(
            broom::augment(econ, se_fit = TRUE)
        )

        expect_equal(
            aug$.fitted,
            unname(predict(ref, type = "link")),
            tolerance = 1e-10
        )

        expect_equal(
            aug$.resid,
            unname(residuals(ref, type = "deviance")),
            tolerance = 1e-10
        )

        expect_equal(
            aug[["(weights)"]],
            unname(weights(ref)),
            tolerance = 1e-10
        )

        expect_equal(
            anyDuplicated(names(aug)),
            0L
        )

        expect_equal(
            nrow(aug),
            nobs(econ)
        )
    }
})





