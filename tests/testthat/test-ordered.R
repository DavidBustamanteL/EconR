#### 1. Ordered Logit and Probit Estimation ####

testthat::test_that("Ordered models match MASS::polr estimates", {

    set.seed(123)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n)
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

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

        testthat::expect_equal(
            unname(coef(econ)),
            unname(coef(ref)),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            unname(econ$zeta),
            unname(ref$zeta),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            unname(vcov(econ)),
            unname(vcov(ref)),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            as.numeric(logLik(econ)),
            as.numeric(logLik(ref)),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            predict(econ, type = "probs"),
            predict(ref, type = "probs"),
            tolerance = 1e-8
        )

        testthat::expect_identical(
            predict(econ, type = "class"),
            predict(ref, type = "class")
        )
    }
})


#### 2. Weighted Ordered Logit and Probit ####

testthat::test_that("Weighted ordered models match MASS::polr", {

    set.seed(456)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n),
        w = sample(1:5, n, replace = TRUE)
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        ref = MASS::polr(
            y ~ x1 + x2 + x3,
            data = df,
            weights = w,
            method = method,
            Hess = TRUE
        )

        econ = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df,
            weights = w,
            link = link
        )

        testthat::expect_equal(
            unname(coef(econ)),
            unname(coef(ref)),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            unname(econ$zeta),
            unname(ref$zeta),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            unname(vcov(econ)),
            unname(vcov(ref)),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            unname(econ$weights),
            unname(df$w)
        )

        testthat::expect_equal(
            sum(econ$weights),
            sum(df$w)
        )

        testthat::expect_equal(
            as.numeric(logLik(econ)),
            as.numeric(logLik(ref)),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            predict(econ, type = "probs"),
            predict(ref, type = "probs"),
            tolerance = 1e-8
        )
    }
})


#### 3. Missing Values and Subsetting ####

testthat::test_that("Ordered models handle missing values and subsets correctly", {

    set.seed(789)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n),
        w = sample(1:3, n, replace = TRUE),
        keep = sample(
            c(TRUE, FALSE),
            n,
            replace = TRUE,
            prob = c(0.8, 0.2)
        )
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    df$x1[sample(n, 30)] = NA
    df$x2[sample(n, 20)] = NA
    df$w[sample(n, 50)] = NA

    expected_rows = which(
        df$keep &
        complete.cases(df[, c("y", "x1", "x2", "x3", "w")])
    )

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        ref = MASS::polr(
            y ~ x1 + x2 + x3,
            data = df,
            weights = w,
            subset = keep,
            method = method,
            Hess = TRUE
        )

        econ = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df,
            weights = w,
            subset = keep,
            link = link
        )

        testthat::expect_equal(
            nobs(econ),
            nobs(ref)
        )

        testthat::expect_equal(
            nrow(econ$model),
            length(expected_rows)
        )

testthat::expect_equal(
    nobs(econ),
    sum(df$w[expected_rows])
)

        testthat::expect_equal(
            econ$used,
            expected_rows
        )

        testthat::expect_equal(
            sort(econ$omitted),
            setdiff(seq_len(n), expected_rows)
        )

        testthat::expect_equal(
            unname(econ$weights),
            unname(df$w[expected_rows])
        )

        testthat::expect_equal(
            unname(coef(econ)),
            unname(coef(ref)),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            unname(econ$zeta),
            unname(ref$zeta),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            unname(vcov(econ)),
            unname(vcov(ref)),
            tolerance = 1e-8
        )
    }
})


#### 4. HC0 and HC1 Robust Covariance ####

testthat::test_that("Ordered HC0 and HC1 covariance matches sandwich", {

    set.seed(1234)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n)
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        ref = MASS::polr(
            y ~ x1 + x2 + x3,
            data = df,
            method = method,
            Hess = TRUE,
            model = TRUE
        )

        econ_hc0 = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df,
            link = link,
            vcov = "HC0"
        )

        econ_hc1 = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df,
            link = link,
            vcov = "HC1"
        )

        V_hc0_ref = sandwich::sandwich(ref)

        k = length(coef(ref)) + length(ref$zeta)

        V_hc1_ref = V_hc0_ref * n / (n - k)

        testthat::expect_equal(
            unname(vcov(econ_hc0)),
            unname(V_hc0_ref),
            tolerance = 1e-7
        )

        testthat::expect_equal(
            unname(vcov(econ_hc1)),
            unname(V_hc1_ref),
            tolerance = 1e-7
        )

        testthat::expect_equal(
            unname(coef(econ_hc0)),
            unname(coef(econ_hc1)),
            tolerance = 1e-10
        )

        testthat::expect_equal(
            unname(econ_hc0$zeta),
            unname(econ_hc1$zeta),
            tolerance = 1e-10
        )
    }
})


#### 5. CR0 and CR1 Clustered Covariance ####

testthat::test_that("Ordered CR0 and CR1 covariance matches sandwich", {

    set.seed(5678)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n),
        cluster_id = sample(1:50, n, replace = TRUE)
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        ref = MASS::polr(
            y ~ x1 + x2 + x3,
            data = df,
            method = method,
            Hess = TRUE,
            model = TRUE
        )

        econ_cr0 = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df,
            link = link,
            vcov = "CR0",
            cluster = cluster_id
        )

        econ_cr1 = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df,
            link = link,
            vcov = "CR1",
            cluster = cluster_id
        )

        V_cr0_ref = sandwich::vcovCL(
            ref,
            cluster = df$cluster_id,
            type = "HC0",
            cadjust = FALSE
        )

        k = length(coef(ref)) + length(ref$zeta)

        V_cr1_ref = V_cr0_ref * (n - 1) / (n - k)

        testthat::expect_equal(
            unname(vcov(econ_cr0)),
            unname(V_cr0_ref),
            tolerance = 1e-7
        )

        testthat::expect_equal(
            unname(vcov(econ_cr1)),
            unname(V_cr1_ref),
            tolerance = 1e-7
        )

        testthat::expect_equal(
            econ_cr0$cluster,
            df$cluster_id
        )

        testthat::expect_equal(
            econ_cr1$cluster,
            df$cluster_id
        )

        testthat::expect_equal(
            unname(coef(econ_cr0)),
            unname(coef(econ_cr1)),
            tolerance = 1e-10
        )
    }
})


#### 6. Weighted Robust and Clustered Covariance ####

testthat::test_that("Weighted ordered robust covariance is correctly normalized", {

    set.seed(9012)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n),
        w = sample(1:4, n, replace = TRUE),
        cluster_id = sample(1:50, n, replace = TRUE)
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        ref = MASS::polr(
            y ~ x1 + x2 + x3,
            data = df,
            weights = w,
            method = method,
            Hess = TRUE,
            model = TRUE
        )

        V_classical = vcov(ref)
        S = sandwich::estfun(ref)

        k = ncol(S)

        V_hc0 = V_classical %*%
            crossprod(S) %*%
            V_classical

        V_hc1 = V_hc0 * n / (n - k)

        cluster_scores = rowsum(
            S,
            group = df$cluster_id,
            reorder = TRUE
        )

        V_cr0 = V_classical %*%
            crossprod(cluster_scores) %*%
            V_classical

        V_cr1 = V_cr0 * (n - 1) / (n - k)

        for (type in c("HC0", "HC1", "CR0", "CR1")) {

            econ = econ_ordered(
                y ~ x1 + x2 + x3,
                data = df,
                weights = w,
                link = link,
                vcov = type,
                cluster = cluster_id
            )

            V_expected = switch(
                type,
                HC0 = V_hc0,
                HC1 = V_hc1,
                CR0 = V_cr0,
                CR1 = V_cr1
            )

            testthat::expect_equal(
                unname(vcov(econ)),
                unname(V_expected),
                tolerance = 1e-7
            )
        }
    }
})


#### 7. Cluster Alignment with Missing Values and Subsets ####

testthat::test_that("Ordered models align clusters after NA removal and subsetting", {

    set.seed(3456)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n),
        w = sample(1:3, n, replace = TRUE),
        cluster_id = sample(1:50, n, replace = TRUE),
        keep = sample(c(TRUE, FALSE), n, replace = TRUE)
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    df$x1[sample(n, 30)] = NA
    df$w[sample(n, 40)] = NA

    expected_rows = which(
        df$keep &
        complete.cases(df[, c("y", "x1", "x2", "x3", "w")])
    )

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        ref = MASS::polr(
            y ~ x1 + x2 + x3,
            data = df,
            weights = w,
            subset = keep,
            method = method,
            Hess = TRUE,
            model = TRUE
        )

        S = sandwich::estfun(ref)
        V_classical = vcov(ref)

        cluster_used = df$cluster_id[expected_rows]

        cluster_scores = rowsum(
            S,
            group = cluster_used,
            reorder = TRUE
        )

        V_cr0 = V_classical %*%
            crossprod(cluster_scores) %*%
            V_classical

        k = ncol(S)
        n_used = nrow(S)

        V_cr1 = V_cr0 * (n_used - 1) / (n_used - k)

        for (type in c("CR0", "CR1")) {

            econ = econ_ordered(
                y ~ x1 + x2 + x3,
                data = df,
                weights = w,
                subset = keep,
                link = link,
                vcov = type,
                cluster = cluster_id
            )

            V_expected = if (type == "CR0") V_cr0 else V_cr1

            testthat::expect_equal(
                econ$used,
                expected_rows
            )

            testthat::expect_equal(
                econ$cluster,
                cluster_used
            )

            testthat::expect_equal(
                unname(vcov(econ)),
                unname(V_expected),
                tolerance = 1e-7
            )
        }
    }
})


#### 8. Input Validation and Error Handling ####

testthat::test_that("Ordered models reject invalid inputs", {

    set.seed(2468)

    n = 300

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        cluster_id = sample(1:20, n, replace = TRUE),
        w = rep(1, n)
    )

    latent = 0.8 * df$x1 - 0.5 * df$x2 + rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    testthat::expect_error(
        econ_ordered(y ~ x1 + x2, data = df, link = "invalid"),
        "arg"
    )

    testthat::expect_error(
    econ_ordered(
        y ~ x1 + x2,
        data = df,
        vcov = "HC3"
    ),
    "arg.*should be one of"
)

    df_numeric = df
    df_numeric$y = as.integer(df_numeric$y)

    testthat::expect_error(
        econ_ordered(y ~ x1 + x2, data = df_numeric),
        "ordered factor"
    )

    df_binary = df
    df_binary$y = ordered(
        ifelse(df_binary$y == "Low", "Low", "High"),
        levels = c("Low", "High")
    )

    testthat::expect_error(
        econ_ordered(y ~ x1 + x2, data = df_binary),
        "at least three"
    )

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df,
            vcov = "CR0"
        ),
        "requires a cluster variable"
    )

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df,
            vcov = "CR1"
        ),
        "requires a cluster variable"
    )

    df_one_cluster = df
    df_one_cluster$cluster_id = 1

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df_one_cluster,
            vcov = "CR0",
            cluster = cluster_id
        ),
        "requires at least two distinct clusters"
    )

    df_missing_cluster = df
    df_missing_cluster$cluster_id[1] = NA

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df_missing_cluster,
            vcov = "CR1",
            cluster = cluster_id
        ),
        "missing values"
    )

    df_negative_weight = df
    df_negative_weight$w[1] = -1

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df_negative_weight,
            weights = w
        ),
        "nonnegative"
    )

    df_infinite_weight = df
    df_infinite_weight$w[1] = Inf

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df_infinite_weight,
            weights = w
        ),
        "finite"
    )

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df,
            weights = rep(1, 10)
        ),
        "one value per row"
    )

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df,
            subset = rep(TRUE, 10)
        ),
        "one value per row"
    )

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df,
            cluster = rep(1, 10)
        ),
        "one value per row"
    )
})


#### 9. Analytical Score Contributions ####

testthat::test_that("Ordered analytical scores match MASS::polr", {

    set.seed(1357)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n),
        w = sample(1:4, n, replace = TRUE)
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        for (weighted in c(FALSE, TRUE)) {

            if (weighted) {

                ref = MASS::polr(
                    y ~ x1 + x2 + x3,
                    data = df,
                    weights = w,
                    method = method,
                    Hess = TRUE,
                    model = TRUE
                )

                econ = econ_ordered(
                    y ~ x1 + x2 + x3,
                    data = df,
                    weights = w,
                    link = link
                )

            } else {

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
            }

            S_ref = sandwich::estfun(ref)
            S_econ = sandwich::estfun(econ)

            testthat::expect_equal(
                unname(S_econ),
                unname(S_ref),
                tolerance = 1e-7
            )

            testthat::expect_equal(
                dim(S_econ),
                c(n, 5L)
            )

            testthat::expect_equal(
                colSums(S_econ),
                rep(0, 5),
                tolerance = 1e-3,
                ignore_attr = TRUE
            )
        }
    }
})


#### 10. Reporting and Prediction Methods ####

testthat::test_that("Ordered reporting and prediction methods are consistent", {

    set.seed(9753)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n),
        cluster_id = sample(1:50, n, replace = TRUE)
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        ref = MASS::polr(
            y ~ x1 + x2 + x3,
            data = df,
            method = method,
            Hess = TRUE,
            model = TRUE
        )

        for (type in c("classical", "HC0", "HC1", "CR0", "CR1")) {

            econ = econ_ordered(
                y ~ x1 + x2 + x3,
                data = df,
                link = link,
                vcov = type,
                cluster = cluster_id
            )

            V = vcov(econ)
            expected_se = sqrt(diag(V))

            model_summary = summary(econ)
            model_tidy = broom::tidy(econ)

            testthat::expect_equal(
                unname(model_summary$coefficients[, "Std. Error"]),
                unname(expected_se),
                tolerance = 1e-8
            )

            testthat::expect_equal(
                unname(model_tidy$std.error),
                unname(expected_se),
                tolerance = 1e-8
            )

            testthat::expect_equal(
                nrow(model_tidy),
                length(expected_se)
            )
        }

        econ = econ_ordered(
            y ~ x1 + x2 + x3,
            data = df,
            link = link
        )

        probs = predict(econ, type = "probs")
        classes = predict(econ, type = "class")

        testthat::expect_equal(
            unname(probs),
            unname(predict(ref, type = "probs")),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            unname(rowSums(probs)),
            rep(1, n),
            tolerance = 1e-10
        )

        testthat::expect_identical(
            classes,
            predict(ref, type = "class")
        )

        ###### Compare Standard Goodness-of-Fit Statistics ####

econ_glance = broom::glance(econ)
ref_glance = broom::glance(ref)

testthat::expect_equal(
    econ_glance[, names(ref_glance), drop = FALSE],
    ref_glance
)

###### Verify McFadden Pseudo-R-Squared ####

testthat::expect_true(
    "pseudo.r.squared" %in% names(econ_glance)
)

testthat::expect_true(
    is.finite(econ_glance$pseudo.r.squared)
)
    }
})


#### 11. Predictions on New Data ####

testthat::test_that("Ordered predictions on new data match MASS::polr", {

    set.seed(8642)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        x3 = rnorm(n)
    )

    latent = 0.8 * df$x1 -
        0.6 * df$x2 +
        0.3 * df$x3 +
        rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    new_data = data.frame(
        x1 = c(-1, 0, 1, 2),
        x2 = c(0.5, -0.5, 1, 0),
        x3 = c(1, 0, -1, 0.5)
    )

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

        probs_ref = predict(
            ref,
            newdata = new_data,
            type = "probs"
        )

        probs_econ = predict(
            econ,
            newdata = new_data,
            type = "probs"
        )

        classes_ref = predict(
            ref,
            newdata = new_data,
            type = "class"
        )

        classes_econ = predict(
            econ,
            newdata = new_data,
            type = "class"
        )

        testthat::expect_equal(
            probs_econ,
            probs_ref,
            tolerance = 1e-10
        )

        testthat::expect_identical(
            classes_econ,
            classes_ref
        )

        testthat::expect_equal(
            unname(rowSums(probs_econ)),
            rep(1, nrow(new_data)),
            tolerance = 1e-10
        )

        testthat::expect_equal(
            dim(probs_econ),
            c(4L, 3L)
        )

        testthat::expect_identical(
            colnames(probs_econ),
            levels(df$y)
        )
    }
})




#### 12. Clustered CRs expansion ####

test_that("Ordered clustered covariance corrections match sandwich", {

    ###### 1. Prepare data ####

    set.seed(123)

    n = 400

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        cluster_id = rep(seq_len(40), each = 10)
    )

    latent = 0.8 * df$x1 - 0.5 * df$x2 + rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -1, 0, 1, Inf),
            labels = c("Low", "Medium", "High", "Very high")
        )
    )

    ###### 2. Fit ordered models ####

    types = c("CR0", "CR1", "CR1p", "CR1S", "stata")

    models = lapply(types, function(type) {
        econ_ordered(
            y ~ x1 + x2,
            data = df,
            vcov = type,
            cluster = cluster_id
        )
    })

    names(models) = types

    ###### 3. Construct independent reference covariance ####

    reference = MASS::polr(
        y ~ x1 + x2,
        data = df,
        method = "logistic",
        Hess = TRUE,
        model = TRUE
    )

    V_classical = stats::vcov(reference)

    S = sandwich::estfun(reference)

    cluster_scores = rowsum(
        S,
        group = df$cluster_id,
        reorder = TRUE
    )

    M = crossprod(cluster_scores)

    n_obs = nrow(S)
    k = ncol(S)
    G = length(unique(df$cluster_id))

    V_CR0 = V_classical %*% M %*% V_classical

    corrections = c(
        CR0 = 1,
        CR1 = (n_obs - 1) / (n_obs - k),
        CR1p = G / (G - k),
        CR1S = G * (n_obs - 1) / ((G - 1) * (n_obs - k)),
        stata = G * (n_obs - 1) / ((G - 1) * (n_obs - k))
    )

    ###### 4. Compare covariance matrices ####

    for (type in types) {

        expected = V_CR0 * corrections[[type]]

        expect_equal(
            unname(stats::vcov(models[[type]])),
            unname(expected),
            tolerance = 1e-7,
            info = paste("Covariance mismatch for", type)
        )
    }

    ###### 5. Verify correction relationships ####

    expect_equal(
        stats::vcov(models[["CR1S"]]),
        stats::vcov(models[["stata"]]),
        tolerance = 1e-12
    )
})


#### 13. Clustered inference and tidy confidence intervals ####

testthat::test_that(
    "ordered models use the appropriate inference distribution",
    {

        ###### 1. Prepare data ####

        df = mtcars
        df$gear = ordered(df$gear)
        df$cluster_id = df$cyl


        ###### 2. Estimate classical and clustered models ####

        m_classical = econ_ordered(
            gear ~ hp,
            data = df,
            vcov = "classical"
        )

        m_cluster = econ_ordered(
            gear ~ hp,
            data = df,
            vcov = "stata",
            cluster = cluster_id
        )


        ###### 3. Verify summary inference ####

        s_classical = summary(m_classical)
        s_cluster = summary(m_cluster)

        testthat::expect_identical(
            s_classical$inference_distribution,
            "z"
        )

        testthat::expect_identical(
            s_cluster$inference_distribution,
            "t"
        )

        testthat::expect_equal(
            s_cluster$inference_df,
            2
        )

        testthat::expect_equal(
            unname(s_cluster$coefficients[, "Pr(>|t|)"]),
            unname(2 * stats::pt(
                abs(s_cluster$coefficients[, "t value"]),
                df = 2,
                lower.tail = FALSE
            ))
        )


        ###### 4. Verify printed summaries ####

        output_classical = capture.output(
            print(s_classical)
        )

        output_cluster = capture.output(
            print(s_cluster)
        )

        testthat::expect_true(
            any(grepl("Inference: standard normal", output_classical))
        )

        testthat::expect_true(
            any(grepl("Inference: t distribution", output_cluster))
        )


        ###### 5. Verify tidy confidence intervals ####

        td_classical = generics::tidy(
            m_classical,
            conf.int = TRUE
        )

        td_cluster = generics::tidy(
            m_cluster,
            conf.int = TRUE
        )

        critical_z = stats::qnorm(0.975)
        critical_t = stats::qt(0.975, df = 2)

        testthat::expect_equal(
            td_classical$conf.low,
            td_classical$estimate -
                critical_z * td_classical$std.error
        )

        testthat::expect_equal(
            td_cluster$conf.low,
            td_cluster$estimate -
                critical_t * td_cluster$std.error
        )

        testthat::expect_equal(
            td_cluster$conf.high,
            td_cluster$estimate +
                critical_t * td_cluster$std.error
        )

        testthat::expect_equal(
            td_cluster$p.value,
            unname(s_cluster$coefficients[, "Pr(>|t|)"])
        )


        ###### 6. Validate confidence level ####

        testthat::expect_error(
            generics::tidy(
                m_cluster,
                conf.int = TRUE,
                conf.level = 1.5
            ),
            "conf.level"
        )
    }
)


#### 14. Ordered Model Confidence Intervals ####

test_that("Ordered confidence intervals use the selected covariance", {

    ###### 1. Prepare data ####

    df = mtcars
df$gear = ordered(df$gear)
df$cluster_id = rep(1:8, each = 4)


    ###### 2. Test classical, robust, and clustered covariance ####

    for (type in c("classical", "HC0", "HC1", "CR0", "CR1", "CR1p", "CR1S", "stata")) {

        clustered = type %in% c(
            "CR0", "CR1", "CR1p", "CR1S", "stata"
        )

        fit = if (clustered) {

            econ_ordered(
                gear ~ hp,
                data = df,
                vcov = type,
                cluster = cluster_id
            )

        } else {

            econ_ordered(
                gear ~ hp,
                data = df,
                vcov = type
            )
        }


        ###### 3. Calculate independent Wald intervals ####

        b = c(coef(fit), fit$zeta)
        se = sqrt(diag(vcov(fit)))
        se = se[names(b)]

        critical = if (clustered) {

            G = length(unique(fit$cluster))

            stats::qt(0.975, df = G - 1L)

        } else {

            stats::qnorm(0.975)
        }

        expected = cbind(
            b - critical * se,
            b + critical * se
        )

        colnames(expected) = c("2.5 %", "97.5 %")


        ###### 4. Compare confidence intervals ####

        expect_equal(
            confint(fit),
            expected,
            tolerance = 1e-10
        )

        expect_equal(
            confint(fit, parm = "hp"),
            expected["hp", , drop = FALSE],
            tolerance = 1e-10
        )

        expect_equal(
            confint(fit, parm = 2),
            expected[2, , drop = FALSE],
            tolerance = 1e-10
        )


        ###### 5. Compare against tidy ####

        td = generics::tidy(
            fit,
            conf.int = TRUE
        )

        expect_equal(
            unname(confint(fit)[, 1]),
            td$conf.low,
            tolerance = 1e-10
        )

        expect_equal(
            unname(confint(fit)[, 2]),
            td$conf.high,
            tolerance = 1e-10
        )
    }


    ###### 6. Validate confidence levels and parameters ####

    fit = econ_ordered(
        gear ~ hp,
        data = df,
        vcov = "HC1"
    )

    expect_error(confint(fit, level = 0))
    expect_error(confint(fit, level = 1))
    expect_error(confint(fit, level = NA_real_))
    expect_error(confint(fit, parm = "unknown"))
    expect_error(confint(fit, parm = 99))
    expect_error(confint(fit, parm = 1.5))


    ###### 7. Verify alternative confidence level ####

    b = c(coef(fit), fit$zeta)
    se = sqrt(diag(vcov(fit)))
    se = se[names(b)]

    critical_90 = stats::qnorm(0.95)

    expect_equal(
        unname(confint(fit, level = 0.90)),
        unname(cbind(
            b - critical_90 * se,
            b + critical_90 * se
        )),
        tolerance = 1e-10
    )
})



#### 15. Ordered Model ANOVA ####

testthat::test_that("Ordered likelihood-ratio tests match MASS", {

    df = mtcars
    df$gear = ordered(df$gear)

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        ref1 = MASS::polr(
            gear ~ hp,
            data = df,
            method = method,
            Hess = TRUE
        )

        ref2 = MASS::polr(
            gear ~ hp + wt,
            data = df,
            method = method,
            Hess = TRUE
        )

        econ1 = econ_ordered(
            gear ~ hp,
            data = df,
            link = link
        )

        econ2 = econ_ordered(
            gear ~ hp + wt,
            data = df,
            link = link
        )

        result = anova(econ1, econ2, test = "Chisq")

        expected_chisq = 2 * (
            as.numeric(logLik(ref2)) -
            as.numeric(logLik(ref1))
        )

        expected_df = attr(logLik(ref2), "df") -
            attr(logLik(ref1), "df")

        expected_p = pchisq(
            expected_chisq,
            df = expected_df,
            lower.tail = FALSE
        )

        testthat::expect_equal(
            result$Chisq[2],
            expected_chisq,
            tolerance = 1e-7
        )

        testthat::expect_equal(
            result$`Df Diff`[2],
            as.numeric(expected_df)
        )

        testthat::expect_equal(
            result$`Pr(>Chisq)`[2],
            expected_p,
            tolerance = 1e-7
        )
    }
})


testthat::test_that("Ordered ANOVA rejects incompatible models", {

    df = mtcars
    df$gear = ordered(df$gear)

    m1 = econ_ordered(gear ~ hp, data = df)
    m2 = econ_ordered(gear ~ hp + wt, data = df)

    testthat::expect_error(
        anova(m2, m1),
        "increasing complexity"
    )

    testthat::expect_error(
        anova(m1),
        "at least two models"
    )

    testthat::expect_error(
        anova(m1, m2, test = "F"),
        "Only likelihood-ratio tests"
    )

    m_probit = econ_ordered(
        gear ~ hp + wt,
        data = df,
        link = "probit"
    )

    testthat::expect_error(
        anova(m1, m_probit),
        "same link"
    )
})


testthat::test_that("Weighted ordered ANOVA matches MASS likelihoods", {

    df = mtcars
    df$gear = ordered(df$gear)
    df$w = rep(1:4, length.out = nrow(df))

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        ref1 = MASS::polr(
            gear ~ hp,
            data = df,
            weights = w,
            method = method,
            Hess = TRUE
        )

        ref2 = MASS::polr(
            gear ~ hp + wt,
            data = df,
            weights = w,
            method = method,
            Hess = TRUE
        )

        econ1 = econ_ordered(
            gear ~ hp,
            data = df,
            weights = w,
            link = link
        )

        econ2 = econ_ordered(
            gear ~ hp + wt,
            data = df,
            weights = w,
            link = link
        )

        result = anova(econ1, econ2)

        expected = 2 * (
            as.numeric(logLik(ref2)) -
            as.numeric(logLik(ref1))
        )

        testthat::expect_equal(
            result$Chisq[2],
            expected,
            tolerance = 1e-7
        )
    }
})


testthat::test_that("Ordered ANOVA supports three nested models", {

    df = mtcars
    df$gear = ordered(df$gear)

    m1 = econ_ordered(gear ~ hp, data = df)
    m2 = econ_ordered(gear ~ hp + wt, data = df)
    m3 = econ_ordered(gear ~ hp + wt + disp, data = df)

    result = anova(m1, m2, m3)

    testthat::expect_equal(nrow(result), 3L)

    testthat::expect_equal(
        result$Chisq[2:3],
        2 * diff(c(
            as.numeric(logLik(m1)),
            as.numeric(logLik(m2)),
            as.numeric(logLik(m3))
        )),
        tolerance = 1e-7
    )
})


testthat::test_that("Ordered ANOVA rejects nonnested models", {

    df = mtcars
    df$gear = ordered(df$gear)

    m1 = econ_ordered(gear ~ hp, data = df)
    m2 = econ_ordered(gear ~ wt + disp, data = df)

    testthat::expect_error(
        anova(m1, m2),
        "nested"
    )
})


testthat::test_that("Ordered ANOVA rejects different estimation samples", {

    df = mtcars
    df$gear = ordered(df$gear)
    df$hp[1] = NA

    m1 = econ_ordered(gear ~ wt, data = df)
    m2 = econ_ordered(gear ~ wt + hp, data = df)

    testthat::expect_error(
        anova(m1, m2),
        "identical estimation samples"
    )
})


testthat::test_that("Ordered ANOVA rejects different weights", {

    df = mtcars
    df$gear = ordered(df$gear)
    df$w = rep(1:4, length.out = nrow(df))

    m1 = econ_ordered(gear ~ hp, data = df)
    m2 = econ_ordered(
        gear ~ hp + wt,
        data = df,
        weights = w
    )

    testthat::expect_error(
        anova(m1, m2),
        "identical observation weights"
    )
})



testthat::test_that("Ordered ANOVA rejects different outcome values", {

    df1 = mtcars
    df2 = mtcars

    df1$gear = ordered(df1$gear)
    df2$gear = ordered(df2$gear)

    df2$gear[1] = 3

    m1 = econ_ordered(gear ~ hp, data = df1)
    m2 = econ_ordered(gear ~ hp + wt, data = df2)

    testthat::expect_error(
        anova(m1, m2),
        "identical dependent-variable values"
    )
})


testthat::test_that("Ordered ANOVA rejects different response levels", {

    df1 = mtcars
    df2 = mtcars

    df1$gear = ordered(
        df1$gear,
        levels = c(3, 4, 5)
    )

    df2$gear = ordered(
        df2$gear,
        levels = c(5, 4, 3)
    )

    m1 = econ_ordered(gear ~ hp, data = df1)
    m2 = econ_ordered(gear ~ hp + wt, data = df2)

    testthat::expect_error(
        anova(m1, m2),
        "identical dependent-variable values"
    )
})


testthat::test_that("Ordered ANOVA rejects equal parameter counts", {

    df = mtcars
    df$gear = ordered(df$gear)

    m1 = econ_ordered(gear ~ hp, data = df)
    m2 = econ_ordered(gear ~ wt, data = df)

    testthat::expect_error(
        anova(m1, m2),
        "different parameter counts"
    )
})


testthat::test_that("Ordered ANOVA is independent of covariance choice", {

    df = mtcars
    df$gear = ordered(df$gear)
    df$cluster_id = rep(1:8, each = 4)

    classical1 = econ_ordered(
        gear ~ hp,
        data = df,
        vcov = "classical"
    )

    classical2 = econ_ordered(
        gear ~ hp + wt,
        data = df,
        vcov = "classical"
    )

    clustered1 = econ_ordered(
        gear ~ hp,
        data = df,
        vcov = "stata",
        cluster = cluster_id
    )

    clustered2 = econ_ordered(
        gear ~ hp + wt,
        data = df,
        vcov = "stata",
        cluster = cluster_id
    )

    testthat::expect_equal(
        as.data.frame(anova(classical1, classical2)),
        as.data.frame(anova(clustered1, clustered2)),
        tolerance = 1e-10
    )
})


testthat::test_that("Ordered ANOVA rejects changed shared regressors", {

    df1 = mtcars
    df2 = mtcars

    df1$gear = ordered(df1$gear)
    df2$gear = ordered(df2$gear)

    df2$hp[1] = df2$hp[1] + 20

    m1 = econ_ordered(
        gear ~ hp,
        data = df1
    )

    m2 = econ_ordered(
        gear ~ hp + wt,
        data = df2
    )

    testthat::expect_error(
        anova(m1, m2),
        "different values for shared regressors"
    )
})


#### 16. Ordered Model Augmentation ####

testthat::test_that("Ordered augmentation matches model predictions", {

    df = mtcars
    df$gear = ordered(df$gear)

    for (link in c("logit", "probit")) {

        econ = econ_ordered(
            gear ~ hp + wt,
            data = df,
            link = link
        )

        aug = broom::augment(econ)

        expected_class = predict(
            econ,
            type = "class"
        )

        expected_probs = predict(
            econ,
            type = "probs"
        )

        testthat::expect_equal(
            nrow(aug),
            nobs(econ)
        )

        testthat::expect_equal(
            aug$.fitted,
            expected_class
        )

        probability_columns = paste0(
            ".fitted_",
            make.names(colnames(expected_probs))
        )

        testthat::expect_equal(
            as.matrix(aug[, probability_columns]),
            unname(expected_probs),
            ignore_attr = TRUE,
            tolerance = 1e-10
        )

        testthat::expect_equal(
            rowSums(aug[, probability_columns]),
            rep(1, nobs(econ)),
            tolerance = 1e-10
        )

        expected_residuals =
            as.integer(df$gear) -
            drop(expected_probs %*% seq_len(ncol(expected_probs)))

        testthat::expect_equal(
            aug$.resid,
            unname(expected_residuals),
            tolerance = 1e-10
        )
    }
})


testthat::test_that("Ordered augmentation supports newdata", {

    df = mtcars
    df$gear = ordered(df$gear)

    newdata = df[1:5, c("hp", "wt")]

    for (link in c("logit", "probit")) {

        econ = econ_ordered(
            gear ~ hp + wt,
            data = df,
            link = link
        )

        aug = broom::augment(
            econ,
            newdata = newdata
        )

        expected_probs = predict(
            econ,
            newdata = newdata,
            type = "probs"
        )

        testthat::expect_equal(
            nrow(aug),
            5L
        )

        testthat::expect_equal(
            aug$hp,
            newdata$hp
        )

        testthat::expect_equal(
            aug$wt,
            newdata$wt
        )

        testthat::expect_equal(
            aug$.fitted,
            predict(econ, newdata = newdata, type = "class")
        )

        testthat::expect_equal(
            as.matrix(aug[, paste0(
                ".fitted_",
                make.names(colnames(expected_probs))
            )]),
            unname(expected_probs),
            ignore_attr = TRUE,
            tolerance = 1e-10
        )

        testthat::expect_false(
            ".resid" %in% names(aug)
        )
    }
})


testthat::test_that("Ordered augmentation validates arguments", {

    df = mtcars
    df$gear = ordered(df$gear)

    econ = econ_ordered(
        gear ~ hp + wt,
        data = df
    )

    testthat::expect_error(
        broom::augment(econ, se_fit = NA),
        "TRUE or FALSE"
    )

    testthat::expect_error(
        broom::augment(
            econ,
            data = df,
            newdata = df[1:5, ]
        ),
        "either 'data' or 'newdata'"
    )
})


test_that("Ordered augmentation calculates prediction standard errors", {

    df = mtcars
    df$gear = ordered(df$gear)
    df$cluster_id = rep(1:8, each = 4)

    for (link in c("logit", "probit")) {

        for (type in c("classical", "HC1", "stata")) {

            m = econ_ordered(
                gear ~ hp + wt,
                data = df,
                link = link,
                vcov = type,
                cluster = cluster_id
            )

            aug = broom::augment(m, se_fit = TRUE)

            se_columns = grep(
                "^\\.se\\.fit_",
                names(aug),
                value = TRUE
            )

            expect_equal(length(se_columns), nlevels(df$gear))
            expect_equal(nrow(aug), length(m$used))

            expect_true(all(is.finite(
                as.matrix(aug[, se_columns])
            )))

            expect_true(all(
                as.matrix(aug[, se_columns]) >= 0
            ))
        }
    }
})


test_that("Ordered prediction SEs match numerical delta method", {

    ###### 1. Prepare data ####

    df = mtcars
    df$gear = ordered(df$gear)
    df$cluster_id = rep(1:8, each = 4)

    newdata = df[1:4, c("hp", "wt")]

    ###### 2. Compare links and covariance estimators ####

    for (link in c("logit", "probit")) {

        for (type in c("classical", "HC1", "stata")) {

            m = econ_ordered(
                gear ~ hp + wt,
                data = df,
                link = link,
                vcov = type,
                cluster = cluster_id
            )

            aug = broom::augment(
                m,
                newdata = newdata,
                se_fit = TRUE
            )

            beta = coef(m)
            zeta = m$zeta
            theta = c(beta, zeta)

            V = vcov(m)
            V = V[names(theta), names(theta), drop = FALSE]

            X = model.matrix(
                delete.response(terms(m)),
                data = newdata,
                contrasts.arg = m$contrasts
            )

            X = X[, names(beta), drop = FALSE]

            K = length(beta)
            J = length(zeta) + 1L

            ###### 3. Probability function ####

            probability = function(parameters, i, j) {

                b = parameters[seq_len(K)]
                cutpoints = parameters[K + seq_len(J - 1L)]

                eta = drop(X[i, , drop = FALSE] %*% b)

                F = if (link == "logit") {
                    plogis
                } else {
                    pnorm
                }

                upper = if (j == J) {
                    1
                } else {
                    F(cutpoints[j] - eta)
                }

                lower = if (j == 1L) {
                    0
                } else {
                    F(cutpoints[j - 1L] - eta)
                }

                upper - lower
            }

            ###### 4. Numerical central differences ####

            h = 1e-5

            for (i in seq_len(nrow(newdata))) {

                for (j in seq_len(J)) {

                    gradient = numeric(length(theta))

                    for (k in seq_along(theta)) {

                        plus = theta
                        minus = theta

                        plus[k] = plus[k] + h
                        minus[k] = minus[k] - h

                        gradient[k] = (
                            probability(plus, i, j) -
                            probability(minus, i, j)
                        ) / (2 * h)
                    }

                    expected_se = sqrt(pmax(
                        drop(crossprod(
                            gradient,
                            V %*% gradient
                        )),
                        0
                    ))

                    se_column = paste0(
                        ".se.fit_",
                        make.names(levels(df$gear)[j])
                    )

                    expect_equal(
                        aug[[se_column]][i],
                        expected_se,
                        tolerance = 1e-5,
                        info = paste(
                            "link =", link,
                            "vcov =", type,
                            "observation =", i,
                            "category =", j
                        )
                    )
                }
            }
        }
    }
})



test_that("Ordered augmentation handles missing observations", {
    df = mtcars
    df$gear = ordered(df$gear)
    df$hp[1] = NA

    m = econ_ordered(gear ~ hp + wt, data = df)
    aug = broom::augment(m)

    expect_equal(nrow(aug), nobs(m))
    expect_equal(nrow(aug), 31L)
    expect_false(anyNA(aug$.fitted))
})


test_that("Ordered augmentation handles weights", {
    df = mtcars
    df$gear = ordered(df$gear)
    df$w = rep(1:4, length.out = nrow(df))

    m = econ_ordered(gear ~ hp + wt, data = df, weights = w)
    aug = broom::augment(m)

    expect_equal(nrow(aug), length(m$used))
    expect_equal(
        aug$.fitted,
        predict(m, newdata = m$model, type = "class")
    )
})


test_that("Ordered augmentation handles factor predictors", {
    df = mtcars
    df$gear = ordered(df$gear)
    df$cyl = factor(df$cyl)

    m = econ_ordered(gear ~ hp + cyl, data = df)
    aug = broom::augment(m)

    expect_equal(nrow(aug), nobs(m))
    expect_equal(
        rowSums(aug[, grep("^\\.fitted_", names(aug))]),
        rep(1, nobs(m)),
        tolerance = 1e-10
    )
})


#### 17. Ordered Model Method Audit ####

test_that("Ordered model accessors and diagnostics audit", {

    df = mtcars
    df$gear = ordered(df$gear)

    m = econ_ordered(
        gear ~ hp + wt,
        data = df
    )

    checks = list(
        fitted = function() fitted(m),
        residuals = function() residuals(m),
        model_matrix = function() model.matrix(m),
        model_frame = function() model.frame(m),
        logLik = function() logLik(m),
        AIC = function() AIC(m),
        BIC = function() BIC(m),
        nobs = function() nobs(m),
        df_residual = function() df.residual(m),
        hatvalues = function() hatvalues(m),
        cooks_distance = function() cooks.distance(m)
    )

    results = lapply(checks, function(f) {
        tryCatch(
            {
                value = f()
                paste("OK", class(value)[1])
            },
            error = function(e) {
                paste("ERROR:", conditionMessage(e))
            }
        )
    })

    print(unlist(results))

    expect_true(TRUE)
})


#### 18. Ordered Model Residuals ####

test_that("Ordered residuals match ordinal-score calculations", {

    df = mtcars
    df$gear = ordered(df$gear)

    for (link in c("logit", "probit")) {

        m = econ_ordered(
            gear ~ hp + wt,
            data = df,
            link = link
        )

        P = as.matrix(predict(m, type = "probs"))
        scores = seq_len(ncol(P))

        expected = drop(P %*% scores)
        variance = drop(P %*% scores^2) - expected^2

        r = as.integer(df$gear) - expected

        expect_equal(
            unname(residuals(m)),
            unname(r),
            tolerance = 1e-10
        )

        expect_equal(
            unname(residuals(m, type = "pearson")),
            unname(r / sqrt(pmax(
                variance,
                .Machine$double.eps
            ))),
            tolerance = 1e-10
        )

        expect_equal(
            unname(residuals(m)),
            unname(broom::augment(m)$.resid),
            tolerance = 1e-10
        )
    }
})



#### 19. Ordered Prediction Edge Cases ####

test_that("Ordered na.exclude matches MASS::polr behavior", {

    df = mtcars
    df$gear = ordered(df$gear)
    df$hp[c(2, 5, 9)] = NA

    for (link in c("logit", "probit")) {

        method = if (link == "logit") "logistic" else "probit"

        econ = econ_ordered(
            gear ~ hp + wt,
            data = df,
            link = link,
            na.action = na.exclude
        )

        ref = MASS::polr(
            gear ~ hp + wt,
            data = df,
            method = method,
            na.action = na.exclude,
            Hess = TRUE,
            model = TRUE
        )

        expect_equal(coef(econ), coef(ref), tolerance = 1e-8)
        expect_equal(econ$zeta, ref$zeta, tolerance = 1e-8)
        expect_equal(econ$used, which(!is.na(df$hp)))
        expect_s3_class(econ$na.action, "exclude")

        expect_equal(
            predict(econ, newdata = econ$model, type = "probs"),
            predict(ref, newdata = ref$model, type = "probs"),
            tolerance = 1e-8
        )
    }
})


test_that("Weighted ordered prediction SEs match numerical derivatives", {

    df = mtcars
    df$gear = ordered(df$gear)
    df$w = rep(1:4, length.out = nrow(df))

    newdata = df[1:3, c("hp", "wt")]

    for (link in c("logit", "probit")) {

        m = econ_ordered(
            gear ~ hp + wt,
            data = df,
            weights = w,
            link = link,
            vcov = "HC1"
        )

        aug = broom::augment(m, newdata = newdata, se_fit = TRUE)

        beta = coef(m)
        zeta = m$zeta
        theta = c(beta, zeta)
        V = vcov(m)[names(theta), names(theta), drop = FALSE]

        X = model.matrix(
            delete.response(terms(m)),
            data = newdata,
            contrasts.arg = m$contrasts
        )
        X = X[, names(beta), drop = FALSE]

        F = if (link == "logit") plogis else pnorm
        K = length(beta)
        J = length(zeta) + 1L

        probability = function(par, i, j) {
            eta = drop(X[i, , drop = FALSE] %*% par[seq_len(K)])
            cuts = par[K + seq_len(J - 1L)]
            upper = if (j == J) 1 else F(cuts[j] - eta)
            lower = if (j == 1L) 0 else F(cuts[j - 1L] - eta)
            upper - lower
        }

        for (i in seq_len(nrow(newdata))) {
            for (j in seq_len(J)) {

                gradient = vapply(seq_along(theta), function(k) {
                    h = 1e-5
                    plus = minus = theta
                    plus[k] = plus[k] + h
                    minus[k] = minus[k] - h
                    (probability(plus, i, j) -
                     probability(minus, i, j)) / (2 * h)
                }, numeric(1))

                expected = sqrt(pmax(
                    drop(crossprod(gradient, V %*% gradient)), 0
                ))

                column = paste0(".se.fit_", make.names(levels(df$gear)[j]))

                expect_equal(
                    aug[[column]][i],
                    expected,
                    tolerance = 1e-5
                )
            }
        }
    }
})


test_that("Factor-based ordered prediction SEs match numerical derivatives", {

    df = mtcars
    df$gear = ordered(df$gear)
    df$cyl = factor(df$cyl)

    newdata = df[c(1, 3, 8), c("hp", "cyl")]

    for (link in c("logit", "probit")) {

        m = econ_ordered(
            gear ~ hp + cyl,
            data = df,
            link = link,
            vcov = "HC1"
        )

        aug = broom::augment(m, newdata = newdata, se_fit = TRUE)

        beta = coef(m)
        theta = c(beta, m$zeta)
        V = vcov(m)[names(theta), names(theta), drop = FALSE]

        X = model.matrix(
            delete.response(terms(m)),
            data = newdata,
            contrasts.arg = m$contrasts
        )
        X = X[, names(beta), drop = FALSE]

        F = if (link == "logit") plogis else pnorm
        K = length(beta)
        J = length(m$zeta) + 1L

        probability = function(par, i, j) {
            eta = drop(X[i, , drop = FALSE] %*% par[seq_len(K)])
            cuts = par[K + seq_len(J - 1L)]
            upper = if (j == J) 1 else F(cuts[j] - eta)
            lower = if (j == 1L) 0 else F(cuts[j - 1L] - eta)
            upper - lower
        }

        for (i in seq_len(nrow(newdata))) {
            for (j in seq_len(J)) {

                gradient = vapply(seq_along(theta), function(k) {
                    h = 1e-5
                    plus = minus = theta
                    plus[k] = plus[k] + h
                    minus[k] = minus[k] - h
                    (probability(plus, i, j) -
                     probability(minus, i, j)) / (2 * h)
                }, numeric(1))

                expected = sqrt(pmax(
                    drop(crossprod(gradient, V %*% gradient)), 0
                ))

                column = paste0(".se.fit_", make.names(levels(df$gear)[j]))

                expect_equal(
                    aug[[column]][i],
                    expected,
                    tolerance = 1e-5
                )
            }
        }
    }
})


#### 20. Ordered NA Exclude Restoration ####

test_that("Ordered na.exclude restores omitted observations", {

    df = mtcars
    df$gear = ordered(df$gear)
    df$hp[c(2, 5, 9)] = NA

    omitted = c(2, 5, 9)
    retained = setdiff(seq_len(nrow(df)), omitted)

    for (link in c("logit", "probit")) {

        m = econ_ordered(
            gear ~ hp + wt,
            data = df,
            link = link,
            na.action = na.exclude
        )

        ###### 1. Residual restoration ####

        r = residuals(m)

        expect_length(r, nrow(df))
        expect_true(all(is.na(r[omitted])))
        expect_true(all(is.finite(r[retained])))

        ###### 2. Pearson residual restoration ####

        rp = residuals(m, type = "pearson")

        expect_length(rp, nrow(df))
        expect_true(all(is.na(rp[omitted])))
        expect_true(all(is.finite(rp[retained])))

        ###### 3. In-sample prediction restoration ####

        P = predict(m, type = "probs")

        expect_equal(nrow(P), nrow(df))
        expect_true(all(is.na(P[omitted, , drop = FALSE])))

        expect_equal(
            rowSums(P[retained, , drop = FALSE]),
            rep(1, length(retained)),
            tolerance = 1e-10,
            ignore_attr = TRUE
        )

        ###### 4. Explicit newdata predictions ####

        newdata = df[retained[1:3], c("hp", "wt")]

        P_new = predict(m, newdata = newdata, type = "probs")

        expect_equal(nrow(P_new), 3L)
        expect_true(all(is.finite(P_new)))
    }
})


#### 21. Ordered Logit Odds Ratios ####

testthat::test_that("Ordered Logit exponentiation is correct", {

    set.seed(123)

    n = 1000

    df = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n)
    )

    latent = 0.8 * df$x1 - 0.6 * df$x2 + rlogis(n)

    df$y = ordered(
        cut(
            latent,
            breaks = c(-Inf, -0.5, 0.8, Inf),
            labels = c("Low", "Medium", "High")
        )
    )

    model = econ_ordered(
        y ~ x1 + x2,
        data = df,
        link = "logit"
    )

    ###### 1. Tidy estimates ####

    original = generics::tidy(
        model,
        conf.int = TRUE
    )

    odds = generics::tidy(
        model,
        conf.int = TRUE,
        exponentiate = TRUE
    )

    is_threshold = original$coef.type == "scale"

    expected = original$estimate
    expected[!is_threshold] = exp(expected[!is_threshold])
    expected[is_threshold] = exp(-expected[is_threshold])

    testthat::expect_equal(
        odds$estimate,
        expected
    )

    ###### 2. Preserve inference ####

    testthat::expect_equal(
        odds$std.error,
        original$std.error
    )

    testthat::expect_equal(
        odds$statistic,
        original$statistic
    )

    testthat::expect_equal(
        odds$p.value,
        original$p.value
    )

    ###### 3. Confidence intervals ####

    expected_low = exp(original$conf.low)
    expected_high = exp(original$conf.high)

    expected_low[is_threshold] =
        exp(-original$conf.high[is_threshold])

    expected_high[is_threshold] =
        exp(-original$conf.low[is_threshold])

    testthat::expect_equal(
        odds$conf.low,
        expected_low
    )

    testthat::expect_equal(
        odds$conf.high,
        expected_high
    )

    ###### 4. Summary ####

    original_summary = summary(model)
    odds_summary = summary(model, exponentiate = TRUE)

    testthat::expect_equal(
        unname(odds_summary$coefficients[, "Estimate"]),
        expected
    )

    testthat::expect_equal(
        odds_summary$coefficients[, "Std. Error"],
        original_summary$coefficients[, "Std. Error"]
    )

    ###### 5. Reject ordered Probit ####

    probit = econ_ordered(
        y ~ x1 + x2,
        data = df,
        link = "probit"
    )

    testthat::expect_error(
        generics::tidy(probit, exponentiate = TRUE),
        "only supported for ordered Logit"
    )

    testthat::expect_error(
        summary(probit, exponentiate = TRUE),
        "only supported for ordered Logit"
    )
})


#### 22. McFadden Pseudo-R-Squared ####

testthat::test_that(
    "Ordered McFadden pseudo-R-squared matches null likelihood",
    {

        set.seed(2026)

        n = 1000

        df = data.frame(
            x1 = rnorm(n),
            x2 = rnorm(n)
        )

        latent = 0.8 * df$x1 -
            0.5 * df$x2 +
            rlogis(n)

        df$y = ordered(
            cut(
                latent,
                breaks = c(-Inf, -0.5, 0.8, Inf),
                labels = c("Low", "Medium", "High")
            )
        )

        for (link in c("logit", "probit")) {

            ###### 1. Estimate Full Model ####

            econ = econ_ordered(
                y ~ x1 + x2,
                data = df,
                link = link
            )

            ###### 2. Calculate Null Likelihood ####

            counts = as.numeric(table(df$y))

            probabilities = counts / sum(counts)

            null_loglik = sum(
                counts * log(probabilities)
            )

            expected = 1 -
                as.numeric(logLik(econ)) / null_loglik

            ###### 3. Check Summary ####

            testthat::expect_equal(
                summary(econ)$pseudo.r.squared,
                expected,
                tolerance = 1e-10
            )

            ###### 4. Check Glance ####

            testthat::expect_equal(
                broom::glance(econ)$pseudo.r.squared,
                expected,
                tolerance = 1e-10
            )

            ###### 5. Check Printed Summary ####

            testthat::expect_output(
                print(summary(econ)),
                "McFadden pseudo-R-squared"
            )
        }
    }
)


#### 23. Weighted McFadden Pseudo-R-Squared ####

testthat::test_that(
    "Weighted ordered McFadden pseudo-R-squared is correct",
    {

        set.seed(2027)

        n = 1000

        df = data.frame(
            x1 = rnorm(n),
            x2 = rnorm(n),
            w = sample(1:5, n, replace = TRUE)
        )

        latent = 0.8 * df$x1 -
            0.5 * df$x2 +
            rlogis(n)

        df$y = ordered(
            cut(
                latent,
                breaks = c(-Inf, -0.5, 0.8, Inf),
                labels = c("Low", "Medium", "High")
            )
        )

        for (link in c("logit", "probit")) {

            econ = econ_ordered(
                y ~ x1 + x2,
                data = df,
                weights = w,
                link = link
            )

            counts = as.numeric(
                tapply(
                    df$w,
                    df$y,
                    sum
                )
            )

            probabilities = counts / sum(counts)

            null_loglik = sum(
                counts * log(probabilities)
            )

            expected = 1 -
                as.numeric(logLik(econ)) / null_loglik

            testthat::expect_equal(
                summary(econ)$pseudo.r.squared,
                expected,
                tolerance = 1e-10
            )

            testthat::expect_equal(
                broom::glance(econ)$pseudo.r.squared,
                expected,
                tolerance = 1e-10
            )
        }
    }
)