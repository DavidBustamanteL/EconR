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
        econ_ordered(y ~ x1 + x2, data = df, vcov = "HC3"),
        "arg"
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
        "require a cluster"
    )

    testthat::expect_error(
        econ_ordered(
            y ~ x1 + x2,
            data = df,
            vcov = "CR1"
        ),
        "require a cluster"
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
        "at least two clusters"
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

        testthat::expect_equal(
            broom::glance(econ),
            broom::glance(ref)
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

