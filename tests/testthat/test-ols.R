#### 1. Basic OLS and Model-Specific NA Handling ####

test_that("OLS matches lm with model-specific NA omission", {

    d = data.frame(
        y = c(1, 2, 3, 4, 5),
        x = c(2, NA, 4, 5, 7),
        irrelevant = c(NA, 1, 2, 3, 4)
    )

    a = econ_lm(
        y ~ x,
        data = d
    )

    b = lm(
        y ~ x,
        data = d
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )

    expect_equal(
        nobs(a),
        nobs(b)
    )

    expect_equal(
        residuals(a),
        residuals(b),
        tolerance = 1e-10
    )

    expect_equal(
        fitted(a),
        fitted(b),
        tolerance = 1e-10
    )
})


#### 2. Formula Features and Interactions ####

test_that("interactions match lm", {

    d = mtcars

    a = econ_lm(
        mpg ~ wt * factor(am),
        data = d
    )

    b = lm(
        mpg ~ wt * factor(am),
        data = d
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )
})


#### 3. Weighted Least Squares ####

test_that("WLS matches lm coefficients", {

    d = mtcars
    w = seq_len(nrow(d))

    a = econ_lm(
        mpg ~ wt + hp,
        data = d,
        weights = w
    )

    b = lm(
        mpg ~ wt + hp,
        data = d,
        weights = w
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )

    expect_equal(
        fitted(a),
        fitted(b),
        tolerance = 1e-10
    )

    expect_equal(
        residuals(a),
        residuals(b),
        tolerance = 1e-10
    )
})


#### 4. HC Covariance Estimators ####

test_that("HC0-HC3 match sandwich lm reference", {

    d = mtcars

    for (type in c(
        "HC0",
        "HC1",
        "HC2",
        "HC3"
    )) {

        a = econ_lm(
            mpg ~ wt + hp,
            data = d,
            vcov = type
        )

        b = lm(
            mpg ~ wt + hp,
            data = d
        )

        ref = sandwich::vcovHC(
            b,
            type = type
        )

        expect_equal(
            unname(as.matrix(vcov(a))),
            unname(as.matrix(ref)),
            tolerance = 1e-8,
            ignore_attr = TRUE
        )

        expect_equal(
            unname(as.matrix(sandwich::vcovHC(a, type = type))),
            unname(as.matrix(ref)),
            tolerance = 1e-8,
            ignore_attr = TRUE
        )
    }
})


#### 5. Missing Weights ####

test_that("missing weights are dropped and warning is issued", {

    d = mtcars %>%
        dplyr::mutate(
            test_weight = seq_len(dplyr::n()),
            test_weight = replace(
                test_weight,
                15,
                NA
            )
        )

    expect_warning(
        a <- econ_lm(
            mpg ~ wt + hp,
            data = d,
            weights = test_weight
        ),
        "Some observations have missingness in the weights variable"
    )

    b = lm(
        mpg ~ wt + hp,
        data = d,
        weights = test_weight
    )

    expect_equal(
        nobs(a),
        nobs(b)
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )
})


#### 6. Clustered Covariance Estimators ####

test_that("CR0, CR2, and stata match estimatr", {

    d = mtcars %>%
        dplyr::mutate(
            cluster_id = rep(1:8, each = 4)
        )

    for (type in c(
        "CR0",
        "CR2",
        "stata"
    )) {

        a = econ_lm(
            mpg ~ wt + hp,
            data = d,
            cluster = cluster_id,
            vcov = type
        )

        b = estimatr::lm_robust(
            mpg ~ wt + hp,
            data = d,
            clusters = cluster_id,
            se_type = type
        )

        expect_equal(
            coef(a),
            coef(b),
            tolerance = 1e-10
        )

        expect_equal(
            unname(as.matrix(vcov(a))),
            unname(as.matrix(vcov(b))),
            tolerance = 1e-8,
            ignore_attr = TRUE
        )
    }
})


#### 7. Clustering with Weights and Missing Values ####

test_that("clustered WLS keeps the correct estimation sample", {

    set.seed(123)

    d = mtcars %>%
        dplyr::mutate(
            cluster_id = rep(1:8, each = 4),
            test_weight = runif(dplyr::n()),
            wt = replace(wt, c(3, 10), NA),
            hp = replace(hp, 5, NA),
            test_weight = replace(
                test_weight,
                15,
                NA
            )
        )

    for (type in c(
        "CR0",
        "CR2",
        "stata"
    )) {

        expect_warning(
            a <- econ_lm(
                mpg ~ wt + hp,
                data = d,
                weights = test_weight,
                cluster = cluster_id,
                vcov = type
            ),
            "Some observations have missingness in the weights variable"
        )

        suppressWarnings(
            b <- estimatr::lm_robust(
                mpg ~ wt + hp,
                data = d,
                weights = test_weight,
                clusters = cluster_id,
                se_type = type
            )
        )

        expect_equal(
            nobs(a),
            28L
        )

        expect_equal(
            nobs(a),
            nobs(b)
        )

        expect_equal(
            coef(a),
            coef(b),
            tolerance = 1e-10
        )

        expect_equal(
            unname(as.matrix(vcov(a))),
            unname(as.matrix(vcov(b))),
            tolerance = 1e-8,
            ignore_attr = TRUE
        )
    }
})


#### 8. Rank Deficiency ####

test_that("rank-deficient HC3 models match estimatr", {

    d = mtcars %>%
        dplyr::mutate(
            wt_copy = wt
        )

    a = econ_lm(
        mpg ~ wt + wt_copy + hp,
        data = d,
        vcov = "HC3"
    )

    b = estimatr::lm_robust(
        mpg ~ wt + wt_copy + hp,
        data = d,
        se_type = "HC3"
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )

    expect_equal(
        unname(as.matrix(vcov(a))),
        unname(as.matrix(vcov(b))),
        tolerance = 1e-8,
        ignore_attr = TRUE
    )

    expect_true(
        is.na(coef(a)["wt_copy"])
    )

    expect_true(
        all(is.na(vcov(a)["wt_copy", ]))
    )

    expect_true(
        all(is.na(vcov(a)[, "wt_copy"]))
    )

    expect_equal(
        a$df.residual,
        29L
    )
})


#### 9. Rank Deficiency Across Covariance Types ####

test_that("rank deficiency works across classical and HC covariance types", {

    d = mtcars %>%
        dplyr::mutate(
            wt_copy = wt
        )

    for (type in c(
        "classical",
        "HC0",
        "HC1",
        "HC2",
        "HC3"
    )) {

        a = econ_lm(
            mpg ~ wt + wt_copy + hp,
            data = d,
            vcov = type
        )

        expect_true(
            is.na(coef(a)["wt_copy"])
        )

        expect_true(
            all(is.na(vcov(a)["wt_copy", ]))
        )

        expect_true(
            all(is.na(vcov(a)[, "wt_copy"]))
        )

        expect_equal(
            a$rank,
            3L
        )

        expect_equal(
            a$df.residual,
            29L
        )
    }
})


#### 10. Rank Deficiency with Clustering ####

test_that("rank-deficient clustered models match estimatr", {

    d = mtcars %>%
        dplyr::mutate(
            wt_copy = wt,
            cluster_id = rep(1:8, each = 4)
        )

    for (type in c(
        "CR0",
        "CR2",
        "stata"
    )) {

        a = econ_lm(
            mpg ~ wt + wt_copy + hp,
            data = d,
            cluster = cluster_id,
            vcov = type
        )

        b = estimatr::lm_robust(
            mpg ~ wt + wt_copy + hp,
            data = d,
            clusters = cluster_id,
            se_type = type
        )

        expect_equal(
            coef(a),
            coef(b),
            tolerance = 1e-10
        )

        expect_equal(
            unname(as.matrix(vcov(a))),
            unname(as.matrix(vcov(b))),
            tolerance = 1e-8,
            ignore_attr = TRUE
        )

        expect_true(
            is.na(coef(a)["wt_copy"])
        )
    }
})


#### 11. Hat Values ####

test_that("hatvalues match lm", {

    d = mtcars

    a = econ_lm(
        mpg ~ wt + hp,
        data = d,
        vcov = "HC3"
    )

    b = lm(
        mpg ~ wt + hp,
        data = d
    )

    expect_equal(
        unname(hatvalues(a)),
        unname(hatvalues(b)),
        tolerance = 1e-10
    )
})


#### 12. Standard Model Methods ####

test_that("standard model methods match lm", {

    d = mtcars

    a = econ_lm(
        mpg ~ wt + hp,
        data = d
    )

    b = lm(
        mpg ~ wt + hp,
        data = d
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )

    expect_equal(
        fitted(a),
        fitted(b),
        tolerance = 1e-10
    )

    expect_equal(
        residuals(a),
        residuals(b),
        tolerance = 1e-10
    )

    expect_equal(
        nobs(a),
        nobs(b)
    )

    newdata = data.frame(
        wt = c(2.5, 3.5),
        hp = c(100, 150)
    )

    expect_equal(
        predict(a, newdata = newdata),
        predict(b, newdata = newdata)
    )
})


#### 13. Degenerate Model Tests ####

test_that("constant outcome produces perfect-fit warning", {

    df_constant = mtcars %>%
        mutate(
            y_constant = 5
        )

    expect_warning(
        econ_lm(
            y_constant ~ wt + hp,
            data = df_constant,
            vcov = "HC1"
        ),
        "Essentially perfect fit: inference may be unreliable."
    )

    m = suppressWarnings(
        econ_lm(
            y_constant ~ wt + hp,
            data = df_constant,
            vcov = "HC1"
        )
    )

    expect_s3_class(
        m,
        "econ_lm"
    )

    expect_equal(
        coef(m)[["(Intercept)"]],
        5,
        tolerance = 1e-12
    )

    expect_equal(
        coef(m)[["wt"]],
        0,
        tolerance = 1e-12
    )

    expect_equal(
        coef(m)[["hp"]],
        0,
        tolerance = 1e-12
    )

    expect_true(
        is.nan(glance.econ_lm(m)$r.squared)
    )
})


test_that("zero usable observations produces informative error", {

    df_all_na = mtcars %>%
        mutate(
            bad_predictor = NA_real_
        )

    expect_error(
        econ_lm(
            mpg ~ wt + bad_predictor,
            data = df_all_na,
            vcov = "HC1"
        ),
        "No complete observations available for estimation."
    )
})


test_that("one usable observation is handled", {

    df_one_obs = mtcars %>%
        mutate(
            x_test = NA_real_,
            x_test = replace(
                x_test,
                1,
                1
            )
        )

    expect_warning(
        econ_lm(
            mpg ~ x_test,
            data = df_one_obs,
            vcov = "HC1"
        ),
        "HC1 covariances become"
    )

    m = suppressWarnings(
        econ_lm(
            mpg ~ x_test,
            data = df_one_obs,
            vcov = "HC1"
        )
    )

    expect_s3_class(
        m,
        "econ_lm"
    )

    expect_equal(
        nobs(m),
        1
    )

    expect_equal(
        m$rank,
        1
    )

    expect_equal(
        m$df.residual,
        0
    )

    expect_equal(
        coef(m)[["(Intercept)"]],
        mtcars$mpg[1]
    )

    expect_true(
        is.na(coef(m)[["x_test"]])
    )

    td = tidy.econ_lm(m)

    expect_true(
        all(is.na(td$conf.low))
    )

    expect_true(
        all(is.na(td$conf.high))
    )
})


test_that("saturated full-rank model is handled", {

    df_exact = mtcars %>%
        slice(1:3)

    expect_warning(
        econ_lm(
            mpg ~ wt + hp,
            data = df_exact,
            vcov = "HC1"
        ),
        "HC1 covariances become"
    )

    m = suppressWarnings(
        econ_lm(
            mpg ~ wt + hp,
            data = df_exact,
            vcov = "HC1"
        )
    )

    expect_s3_class(
        m,
        "econ_lm"
    )

    expect_equal(
        nobs(m),
        3
    )

    expect_equal(
        m$rank,
        3
    )

    expect_equal(
        m$df.residual,
        0
    )

    expect_false(
        any(is.na(coef(m)))
    )

    td = tidy.econ_lm(m)

    expect_true(
        all(is.na(td$conf.low))
    )

    expect_true(
        all(is.na(td$conf.high))
    )
})

#### 14. All-Zero Weights ####

test_that("all-zero weights produce informative error", {

    df_zero_weights = mtcars %>%
        mutate(
            test_weight = 0
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = df_zero_weights,
            weights = test_weight,
            vcov = "HC1"
        ),
        "At least one observation must have a positive weight."
    )
})


#### 15. WLS Goodness-of-Fit ####

test_that("WLS goodness-of-fit matches lm with positive unequal weights", {

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

    g_econ = glance(m_econ)
    s_lm = summary(m_lm)

    expect_equal(
        g_econ$r.squared,
        unname(s_lm$r.squared),
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$adj.r.squared,
        unname(s_lm$adj.r.squared),
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$aic,
        unname(AIC(m_lm)),
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$bic,
        unname(BIC(m_lm)),
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$nobs,
        nobs(m_lm)
    )
})


test_that("WLS goodness-of-fit matches lm with zero weights", {

    df_wls_zero = mtcars %>%
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
        mpg ~ wt + hp,
        data = df_wls_zero,
        weights = test_weight,
        vcov = "HC1"
    )

    m_lm = lm(
        mpg ~ wt + hp,
        data = df_wls_zero,
        weights = test_weight
    )

    g_econ = glance(m_econ)
    s_lm = summary(m_lm)

    expect_equal(
        g_econ$r.squared,
        unname(s_lm$r.squared),
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$adj.r.squared,
        unname(s_lm$adj.r.squared),
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$aic,
        unname(AIC(m_lm)),
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$bic,
        unname(BIC(m_lm)),
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$nobs,
        nobs(m_lm)
    )

    expect_equal(
        m_econ$df.residual,
        m_lm$df.residual
    )
})


test_that("zero-weight observations use effective sample size in sandwich bread", {

    df_zero = mtcars %>%
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
        mpg ~ wt + hp,
        data = df_zero,
        weights = test_weight,
        vcov = "HC1"
    )

    m_lm = lm(
        mpg ~ wt + hp,
        data = df_zero,
        weights = test_weight
    )

    expect_equal(
        nobs(m_econ),
        nobs(m_lm)
    )

    expect_equal(
        m_econ$df.residual,
        m_lm$df.residual
    )

    expect_equal(
        unname(
            sandwich::bread(m_econ)
        ),
        unname(
            sandwich::bread(m_lm)
        ),
        tolerance = 1e-12
    )

    expect_equal(
        unname(vcov(m_econ)),
        unname(
            sandwich::vcovHC(
                m_lm,
                type = "HC1"
            )
        ),
        tolerance = 1e-12
    )
})


test_that("zero weights use effective sample size for WLS reporting", {

    df_wls_zero = mtcars %>%
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
        mpg ~ wt + hp,
        data = df_wls_zero,
        weights = test_weight,
        vcov = "HC1"
    )

    m_lm = lm(
        mpg ~ wt + hp,
        data = df_wls_zero,
        weights = test_weight
    )

    g_econ = generics::glance(m_econ)
    g_lm = broom::glance(m_lm)

    expect_equal(
        nobs(m_econ),
        nobs(m_lm)
    )

    expect_equal(
        m_econ$df.residual,
        df.residual(m_lm)
    )

    expect_equal(
        g_econ$r.squared,
        g_lm$r.squared,
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$adj.r.squared,
        g_lm$adj.r.squared,
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$aic,
        g_lm$AIC,
        tolerance = 1e-12
    )

    expect_equal(
        g_econ$bic,
        g_lm$BIC,
        tolerance = 1e-12
    )
})


#### 16. NA Exclude Semantics ####

test_that("na.exclude restores omitted observations in model outputs", {

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

    expect_equal(
        nobs(m_econ),
        nobs(m_lm)
    )

    expect_equal(
        fitted(m_econ),
        fitted(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        residuals(m_econ),
        residuals(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        predict(m_econ),
        predict(m_lm),
        tolerance = 1e-12
    )
})


test_that("na.omit does not restore omitted observations", {

    df_na_omit = mtcars %>%
        mutate(
            wt = replace(
                wt,
                c(3, 8, 15),
                NA
            )
        )

    m_econ = econ_lm(
        mpg ~ wt + hp,
        data = df_na_omit,
        na.action = na.omit,
        vcov = "HC1"
    )

    m_lm = lm(
        mpg ~ wt + hp,
        data = df_na_omit,
        na.action = na.omit
    )

    expect_equal(
        fitted(m_econ),
        fitted(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        residuals(m_econ),
        residuals(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        predict(m_econ),
        predict(m_lm),
        tolerance = 1e-12
    )
})


#### 17. Offset Terms ####

test_that("offset models match lm", {

    m_econ = econ_lm(
        mpg ~ wt + hp + offset(disp / 100),
        data = mtcars,
        vcov = "HC1"
    )

    m_lm = lm(
        mpg ~ wt + hp + offset(disp / 100),
        data = mtcars
    )

    expect_equal(
        coef(m_econ),
        coef(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        fitted(m_econ),
        fitted(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        residuals(m_econ),
        residuals(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        unname(vcov(m_econ)),
        unname(
            sandwich::vcovHC(
                m_lm,
                type = "HC1"
            )
        ),
        tolerance = 1e-12
    )

    expect_equal(
        unname(
            sandwich::estfun(m_econ)
        ),
        unname(
            sandwich::estfun(m_lm)
        ),
        tolerance = 1e-12
    )
})


test_that("offset prediction evaluates offset in new data", {

    m_econ = econ_lm(
        mpg ~ wt + hp + offset(disp / 100),
        data = mtcars,
        vcov = "HC1"
    )

    m_lm = lm(
        mpg ~ wt + hp + offset(disp / 100),
        data = mtcars
    )

    new_offset = mtcars %>%
        slice(1:5) %>%
        mutate(
            wt = wt * 1.05,
            disp = disp * 1.10
        )

    expect_equal(
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
})


#### 18. Custom Contrasts ####

test_that("custom sum contrasts match lm", {

    df_contrast = mtcars %>%
        mutate(
            cyl_factor = factor(cyl)
        )

    contrasts(df_contrast$cyl_factor) = contr.sum(3)

    m_econ = econ_lm(
        mpg ~ wt + cyl_factor,
        data = df_contrast,
        vcov = "HC1"
    )

    m_lm = lm(
        mpg ~ wt + cyl_factor,
        data = df_contrast
    )

    expect_equal(
        coef(m_econ),
        coef(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        model.matrix(m_econ),
        model.matrix(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        fitted(m_econ),
        fitted(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        residuals(m_econ),
        residuals(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        unname(vcov(m_econ)),
        unname(
            sandwich::vcovHC(
                m_lm,
                type = "HC1"
            )
        ),
        tolerance = 1e-12
    )

    newdata = data.frame(
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

    expect_equal(
        predict(
            m_econ,
            newdata = newdata
        ),
        predict(
            m_lm,
            newdata = newdata
        ),
        tolerance = 1e-12
    )
})


test_that("ordered-factor polynomial contrasts match lm", {

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

    expect_equal(
        coef(m_econ),
        coef(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        model.matrix(m_econ),
        model.matrix(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        unname(vcov(m_econ)),
        unname(
            sandwich::vcovHC(
                m_lm,
                type = "HC1"
            )
        ),
        tolerance = 1e-12
    )

    newdata = data.frame(
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

    expect_equal(
        predict(
            m_econ,
            newdata = newdata
        ),
        predict(
            m_lm,
            newdata = newdata
        ),
        tolerance = 1e-12
    )
})


#### 19. Weighted Offset Models ####

test_that("weighted offset models match lm", {

    df_offset = mtcars %>%
        mutate(
            test_weight = seq(
                0.5,
                2,
                length.out = n()
            )
        )

    m_econ = econ_lm(
        mpg ~ wt + hp + offset(disp / 100),
        data = df_offset,
        weights = test_weight,
        vcov = "HC1"
    )

    m_lm = lm(
        mpg ~ wt + hp + offset(disp / 100),
        data = df_offset,
        weights = test_weight
    )

    expect_equal(
        coef(m_econ),
        coef(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        fitted(m_econ),
        fitted(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        residuals(m_econ),
        residuals(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        unname(vcov(m_econ)),
        unname(
            sandwich::vcovHC(
                m_lm,
                type = "HC1"
            )
        ),
        tolerance = 1e-12
    )
})


test_that("weighted offset models handle missingness correctly", {

    df_offset = mtcars %>%
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
            data = df_offset,
            weights = test_weight,
            vcov = "HC1"
        )
    )

    m_lm = lm(
        mpg ~ wt + hp + offset(offset_var),
        data = df_offset,
        weights = test_weight
    )

    expect_equal(
        nobs(m_econ),
        nobs(m_lm)
    )

    expect_equal(
        coef(m_econ),
        coef(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        fitted(m_econ),
        fitted(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        residuals(m_econ),
        residuals(m_lm),
        tolerance = 1e-12
    )

    expect_equal(
        unname(vcov(m_econ)),
        unname(
            sandwich::vcovHC(
                m_lm,
                type = "HC1"
            )
        ),
        tolerance = 1e-12
    )
})


#### 20. NA Fail Semantics ####

test_that("na.fail matches lm behavior", {

    #### Missing covariate ####

    df_na_fail = mtcars %>%
        mutate(
            hp = replace(
                hp,
                c(3, 8),
                NA
            )
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = df_na_fail,
            na.action = na.fail
        ),
        "missing values in object"
    )


    #### Missing weights ####

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

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = df_na_fail_weight,
            weights = test_weight,
            na.action = na.fail
        ),
        "missing values in object"
    )


    #### Missing cluster IDs ####

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

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = df_na_fail_cluster,
            cluster = cluster_id,
            vcov = "CR0",
            na.action = na.fail
        ),
        "missing values in object"
    )


    #### Complete data ####

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

    expect_equal(
        coef(m_econ),
        coef(m_lm)
    )

    expect_equal(
        fitted(m_econ),
        fitted(m_lm)
    )

    expect_equal(
        residuals(m_econ),
        residuals(m_lm)
    )

    expect_equal(
        nobs(m_econ),
        nobs(m_lm)
    )
})


#### 21. Subset Edge Cases ####

test_that("subset handling matches expected estimation samples", {

    #### Basic subset ####

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

    expect_equal(
        coef(m_econ),
        coef(m_lm)
    )

    expect_equal(
        fitted(m_econ),
        fitted(m_lm)
    )

    expect_equal(
        residuals(m_econ),
        residuals(m_lm)
    )

    expect_equal(
        nobs(m_econ),
        nobs(m_lm)
    )


    #### Missing covariate outside subset ####

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

    expect_equal(
        coef(m_econ),
        coef(m_lm)
    )

    expect_equal(
        nobs(m_econ),
        7L
    )


    #### Missing covariate inside subset ####

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

    expect_equal(
        coef(m_econ),
        coef(m_lm)
    )

    expect_equal(
        rownames(model.frame(m_econ)),
        rownames(model.frame(m_lm))
    )


    #### Missing weight outside subset ####

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

    expect_silent(
        econ_lm(
            mpg ~ wt + hp,
            data = df_subset_weight,
            subset = cyl == 6,
            weights = test_weight
        )
    )


    #### Missing weight inside subset ####

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

    expect_warning(
        econ_lm(
            mpg ~ wt + hp,
            data = df_subset_weight_inside,
            subset = cyl == 6,
            weights = test_weight
        ),
        "missingness in the weights"
    )


    #### Missing cluster outside subset ####

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

    expect_silent(
        econ_lm(
            mpg ~ wt + hp,
            data = df_subset_cluster,
            subset = cyl == 6,
            cluster = cluster_id,
            vcov = "CR0"
        )
    )


    #### Missing cluster inside subset ####

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

    expect_warning(
        econ_lm(
            mpg ~ wt + hp,
            data = df_subset_cluster_inside,
            subset = cyl == 6,
            cluster = cluster_id,
            vcov = "CR0"
        ),
        "missingness in the cluster"
    )


    #### na.fail outside subset ####

    df_subset_fail = mtcars %>%
        mutate(
            hp = replace(
                hp,
                which(cyl != 6)[1],
                NA
            )
        )

    expect_silent(
        econ_lm(
            mpg ~ wt + hp,
            data = df_subset_fail,
            subset = cyl == 6,
            na.action = na.fail
        )
    )


    #### na.fail inside subset ####

    df_subset_fail_inside = mtcars %>%
        mutate(
            hp = replace(
                hp,
                which(cyl == 6)[1],
                NA
            )
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = df_subset_fail_inside,
            subset = cyl == 6,
            na.action = na.fail
        ),
        "missing values in object"
    )
})


#### 22. Formula and Environment Edge Cases ####

test_that("formula and auxiliary inputs respect their environments", {

    #### External predictor ####

    external_x = mtcars$wt * 2

    m_econ = econ_lm(
        mpg ~ external_x + hp,
        data = mtcars
    )

    m_lm = lm(
        mpg ~ external_x + hp,
        data = mtcars
    )

    expect_equal(
        coef(m_econ),
        coef(m_lm)
    )

    expect_equal(
        model.matrix(m_econ),
        model.matrix(m_lm)
    )


    #### Custom formula function ####

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

    expect_equal(
        coef(m_econ),
        coef(m_lm)
    )

    expect_equal(
        model.matrix(m_econ),
        model.matrix(m_lm)
    )


    #### External weights ####

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

    expect_equal(
        coef(m_econ),
        coef(m_lm)
    )

    expect_equal(
        weights(m_econ),
        weights(m_lm)
    )


    #### External subset ####

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

    expect_equal(
        coef(m_econ),
        coef(m_lm)
    )

    expect_equal(
        rownames(model.frame(m_econ)),
        rownames(model.frame(m_lm))
    )


    #### External cluster ####

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

    expect_equal(
        coef(m_econ),
        coef(m_internal)
    )

    expect_equal(
        vcov(m_econ),
        vcov(m_internal)
    )


    #### Local formula environment ####

    fit_local_models = function() {

        local_multiplier = 2.5

        local_transform = function(x) {
            x * local_multiplier
        }

        list(
            econ = econ_lm(
                mpg ~ local_transform(wt) + hp,
                data = mtcars
            ),
            lm = lm(
                mpg ~ local_transform(wt) + hp,
                data = mtcars
            )
        )
    }

    local_models = fit_local_models()

    expect_equal(
        coef(local_models$econ),
        coef(local_models$lm)
    )

    expect_equal(
        model.matrix(local_models$econ),
        model.matrix(local_models$lm)
    )

    expect_equal(
        fitted(local_models$econ),
        fitted(local_models$lm)
    )

    expect_equal(
        residuals(local_models$econ),
        residuals(local_models$lm)
    )
})


#### 23. Prediction Edge Cases ####

test_that("prediction handles edge cases consistently with lm", {

    #### Missing predictor values ####

    m_econ = econ_lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    m_lm = lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    newdata_na = data.frame(
        wt = c(2.5, NA, 3.5),
        hp = c(100, 120, NA)
    )

    expect_equal(
        predict(
            m_econ,
            newdata = newdata_na
        ),
        predict(
            m_lm,
            newdata = newdata_na
        )
    )


    #### Missing predictor column ####

    newdata_missing = data.frame(
        wt = c(2.5, 3.0, 3.5)
    )

    expect_error(
        predict(
            m_econ,
            newdata = newdata_missing
        ),
        "object 'hp' not found"
    )


    #### Factor levels ####

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
            levels = c(4, 6, 8, 10)
        )
    )

    expect_error(
        predict(
            m_econ_factor,
            newdata = newdata_new_level
        ),
        "factor cyl_factor has new level 10"
    )

    newdata_factor = data.frame(
        wt = c(2.5, 3.0, 3.5),
        cyl_factor = factor(
            c(4, 6, 8),
            levels = c(4, 6, 8)
        )
    )

    expect_equal(
        predict(
            m_econ_factor,
            newdata = newdata_factor
        ),
        predict(
            m_lm_factor,
            newdata = newdata_factor
        )
    )


    #### Transformed terms ####

    m_econ_transform = econ_lm(
        mpg ~ log(wt) + I(hp^2),
        data = mtcars
    )

    m_lm_transform = lm(
        mpg ~ log(wt) + I(hp^2),
        data = mtcars
    )

    newdata_transform = data.frame(
        wt = c(2.5, 3.0, 3.5),
        hp = c(100, 150, 200)
    )

    expect_equal(
        predict(
            m_econ_transform,
            newdata = newdata_transform
        ),
        predict(
            m_lm_transform,
            newdata = newdata_transform
        )
    )


    #### Rank-deficient model ####

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
        wt = c(2.5, 3.0, 3.5),
        wt_copy = c(2.5, 3.0, 3.5),
        hp = c(100, 150, 200)
    )

    expect_equal(
        predict(
            m_econ_rank,
            newdata = newdata_rank
        ),
        suppressWarnings(
            predict(
                m_lm_rank,
                newdata = newdata_rank
            )
        )
    )


    #### Interaction terms ####

    m_econ_interaction = econ_lm(
        mpg ~ wt * hp,
        data = mtcars
    )

    m_lm_interaction = lm(
        mpg ~ wt * hp,
        data = mtcars
    )

    newdata_interaction = data.frame(
        wt = c(2.5, 3.0, 3.5),
        hp = c(100, 150, 200)
    )

    expect_equal(
        predict(
            m_econ_interaction,
            newdata = newdata_interaction
        ),
        predict(
            m_lm_interaction,
            newdata = newdata_interaction
        )
    )


    #### na.omit in newdata ####

    newdata_omit = data.frame(
        wt = c(2.5, NA, 3.5),
        hp = c(100, 120, 200)
    )

    expect_equal(
        predict(
            m_econ,
            newdata = newdata_omit,
            na.action = na.omit
        ),
        predict(
            m_lm,
            newdata = newdata_omit,
            na.action = na.omit
        )
    )
})


#### 24. Degrees of Freedom and Tiny Samples ####

test_that(
    "saturated models with zero residual degrees of freedom match lm",
    {

        df_test = data.frame(
            y = c(2, 5, 9),
            x1 = c(1, 2, 4),
            x2 = c(2, 5, 3)
        )

        m_econ = econ_lm(
            y ~ x1 + x2,
            data = df_test
        )

        m_lm = lm(
            y ~ x1 + x2,
            data = df_test
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            df.residual(m_econ),
            0
        )

        expect_equal(
            fitted(m_econ),
            fitted(m_lm)
        )

        expect_equal(
            residuals(m_econ),
            residuals(m_lm)
        )

        expect_equal(
            summary(m_econ)$coefficients,
            summary(m_lm)$coefficients
        )

        expect_true(
            all(
                is.nan(
                    summary(m_econ)$coefficients[
                        ,
                        "Std. Error"
                    ]
                )
            )
        )
    }
)


test_that(
    "models with one residual degree of freedom match lm inference",
    {

        df_test = data.frame(
            y = c(2, 5, 9, 10),
            x1 = c(1, 2, 4, 5),
            x2 = c(2, 5, 3, 7)
        )

        m_econ = econ_lm(
            y ~ x1 + x2,
            data = df_test
        )

        m_lm = lm(
            y ~ x1 + x2,
            data = df_test
        )

        expect_equal(
            df.residual(m_econ),
            1
        )

        expect_equal(
            summary(m_econ)$coefficients,
            summary(m_lm)$coefficients
        )
    }
)


test_that(
    "single-observation intercept-only models match lm",
    {

        df_test = data.frame(
            y = 5
        )

        m_econ = econ_lm(
            y ~ 1,
            data = df_test
        )

        m_lm = lm(
            y ~ 1,
            data = df_test
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            df.residual(m_econ),
            0
        )

        expect_equal(
            fitted(m_econ),
            fitted(m_lm)
        )

        expect_equal(
            residuals(m_econ),
            residuals(m_lm)
        )

        expect_equal(
            summary(m_econ)$coefficients,
            summary(m_lm)$coefficients
        )
    }
)


test_that(
    "constant outcomes return undefined R squared",
    {

        df_test = mtcars %>%
            dplyr::mutate(
                constant_y = 10
            )

        expect_warning(
            econ_lm(
                constant_y ~ wt + hp,
                data = df_test
            ),
            "Essentially perfect fit"
        )

        m_econ = suppressWarnings(
            econ_lm(
                constant_y ~ wt + hp,
                data = df_test
            )
        )

        g = glance(
            m_econ
        )

        expect_true(
            is.nan(
                g$r.squared
            )
        )

        expect_true(
            is.nan(
                g$adj.r.squared
            )
        )

        expect_equal(
            g$rmse,
            0
        )
    }
)


test_that(
    "constant predictors are handled as aliased terms",
    {

        df_test = mtcars %>%
            dplyr::mutate(
                constant_x = 5
            )

        m_econ = econ_lm(
            mpg ~ constant_x + wt,
            data = df_test
        )

        m_lm = lm(
            mpg ~ constant_x + wt,
            data = df_test
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            m_econ$rank,
            m_lm$rank
        )

        expect_equal(
            df.residual(m_econ),
            df.residual(m_lm)
        )

        expect_equal(
            fitted(m_econ),
            fitted(m_lm)
        )
    }
)


test_that(
    "intercept-only models match lm inference",
    {

        m_econ = econ_lm(
            mpg ~ 1,
            data = mtcars
        )

        m_lm = lm(
            mpg ~ 1,
            data = mtcars
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            df.residual(m_econ),
            df.residual(m_lm)
        )

        expect_equal(
            summary(m_econ)$coefficients,
            summary(m_lm)$coefficients
        )
    }
)


test_that(
    "no-intercept models use correct goodness-of-fit definitions",
    {

        m_econ = econ_lm(
            mpg ~ 0 + wt + hp,
            data = mtcars
        )

        m_lm = lm(
            mpg ~ 0 + wt + hp,
            data = mtcars
        )

        g = glance(
            m_econ
        )

        s = summary(
            m_lm
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            g$r.squared,
            s$r.squared
        )

        expect_equal(
            g$adj.r.squared,
            s$adj.r.squared
        )
    }
)


test_that(
    "tidy returns NaN inference quantities when residual df is zero",
    {

        df_test = data.frame(
            y = c(2, 5, 9),
            x1 = c(1, 2, 4),
            x2 = c(2, 5, 3)
        )

        m_econ = econ_lm(
            y ~ x1 + x2,
            data = df_test
        )

        td = tidy(
            m_econ
        )

        expect_true(
            all(
                is.nan(
                    td$std.error
                )
            )
        )

        expect_true(
            all(
                is.nan(
                    td$statistic
                )
            )
        )

        expect_true(
            all(
                is.nan(
                    td$p.value
                )
            )
        )

        expect_true(
            all(
                is.nan(
                    td$conf.low
                )
            )
        )

        expect_true(
            all(
                is.nan(
                    td$conf.high
                )
            )
        )
    }
)


test_that(
    "weighted models with one residual degree of freedom match lm",
    {

        df_test = data.frame(
            y = c(2, 5, 9, 10),
            x1 = c(1, 2, 4, 5),
            x2 = c(2, 5, 3, 7),
            w = c(1, 2, 3, 4)
        )

        m_econ = econ_lm(
            y ~ x1 + x2,
            data = df_test,
            weights = w
        )

        m_lm = lm(
            y ~ x1 + x2,
            data = df_test,
            weights = w
        )

        expect_equal(
            df.residual(m_econ),
            1
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            summary(m_econ)$coefficients,
            summary(m_lm)$coefficients
        )
    }
)


test_that(
    "weighted saturated models handle zero residual df",
    {

        df_test = data.frame(
            y = c(2, 5, 9),
            x1 = c(1, 2, 4),
            x2 = c(2, 5, 3),
            w = c(1, 2, 3)
        )

        m_econ = econ_lm(
            y ~ x1 + x2,
            data = df_test,
            weights = w
        )

        m_lm = lm(
            y ~ x1 + x2,
            data = df_test,
            weights = w
        )

        expect_equal(
            df.residual(m_econ),
            0
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            summary(m_econ)$coefficients,
            summary(m_lm)$coefficients
        )
    }
)


test_that(
    "zero weights determine effective sample size and residual df",
    {

        df_test = data.frame(
            y = c(2, 5, 9, 10, 14),
            x = c(1, 2, 4, 5, 7),
            w = c(1, 1, 1, 0, 0)
        )

        m_econ = econ_lm(
            y ~ x,
            data = df_test,
            weights = w
        )

        m_lm = lm(
            y ~ x,
            data = df_test,
            weights = w
        )

        expect_equal(
            nobs(m_econ),
            nobs(m_lm)
        )

        expect_equal(
            df.residual(m_econ),
            df.residual(m_lm)
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            summary(m_econ)$coefficients,
            summary(m_lm)$coefficients,
            tolerance = 1e-12
        )
    }
)


#### 25. Numerical and Pathological Inputs ####


test_that(
    "very large predictors remain numerically stable",
    {

        df_test = mtcars %>%
            dplyr::mutate(
                wt_large = wt * 1e10,
                hp_large = hp * 1e8
            )

        m_econ = econ_lm(
            mpg ~ wt_large + hp_large,
            data = df_test
        )

        m_lm = lm(
            mpg ~ wt_large + hp_large,
            data = df_test
        )

        expect_equal(
            m_econ$rank,
            m_lm$rank
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm),
            tolerance = 1e-10
        )

        expect_equal(
            fitted(m_econ),
            fitted(m_lm),
            tolerance = 1e-10
        )

        expect_equal(
            vcov(m_econ),
            vcov(m_lm),
            tolerance = 1e-10
        )
    }
)


test_that(
    "very small predictors remain numerically stable",
    {

        df_test = mtcars %>%
            dplyr::mutate(
                wt_small = wt * 1e-10,
                hp_small = hp * 1e-8
            )

        m_econ = econ_lm(
            mpg ~ wt_small + hp_small,
            data = df_test
        )

        m_lm = lm(
            mpg ~ wt_small + hp_small,
            data = df_test
        )

        expect_equal(
            m_econ$rank,
            m_lm$rank
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm),
            tolerance = 1e-10
        )

        expect_equal(
            fitted(m_econ),
            fitted(m_lm),
            tolerance = 1e-10
        )

        expect_equal(
            vcov(m_econ),
            vcov(m_lm),
            tolerance = 1e-10
        )
    }
)


test_that(
    "extreme outcome scaling remains numerically stable",
    {

        df_test = mtcars %>%
            dplyr::mutate(
                mpg_large = mpg * 1e12,
                mpg_small = mpg * 1e-12
            )

        for (response in c(
            "mpg_large",
            "mpg_small"
        )) {

            formula_test = stats::reformulate(
                c(
                    "wt",
                    "hp"
                ),
                response = response
            )

            m_econ = econ_lm(
                formula_test,
                data = df_test
            )

            m_lm = lm(
                formula_test,
                data = df_test
            )

            expect_equal(
                coef(m_econ),
                coef(m_lm),
                tolerance = 1e-10
            )

            expect_equal(
                fitted(m_econ),
                fitted(m_lm),
                tolerance = 1e-10
            )

            expect_equal(
                residuals(m_econ),
                residuals(m_lm),
                tolerance = 1e-10
            )

            expect_equal(
                vcov(m_econ),
                vcov(m_lm),
                tolerance = 1e-10
            )
        }
    }
)


test_that(
    "highly heterogeneous weights match lm",
    {

        df_test = mtcars %>%
            dplyr::mutate(
                w_extreme = 10^seq(
                    -6,
                    6,
                    length.out = n()
                )
            )

        m_econ = econ_lm(
            mpg ~ wt + hp,
            data = df_test,
            weights = w_extreme
        )

        m_lm = lm(
            mpg ~ wt + hp,
            data = df_test,
            weights = w_extreme
        )

        expect_equal(
            m_econ$rank,
            m_lm$rank
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm),
            tolerance = 1e-10
        )

        expect_equal(
            fitted(m_econ),
            fitted(m_lm),
            tolerance = 1e-10
        )

        expect_equal(
            residuals(m_econ),
            residuals(m_lm),
            tolerance = 1e-10
        )

        expect_equal(
            vcov(m_econ),
            vcov(m_lm),
            tolerance = 1e-10
        )
    }
)


test_that(
    "zero and highly unbalanced weights match lm",
    {

        df_test = mtcars %>%
            dplyr::mutate(
                w_mix = c(
                    rep(
                        0,
                        4
                    ),
                    10^seq(
                        -8,
                        8,
                        length.out = n() - 4
                    )
                )
            )

        m_econ = econ_lm(
            mpg ~ wt + hp,
            data = df_test,
            weights = w_mix
        )

        m_lm = lm(
            mpg ~ wt + hp,
            data = df_test,
            weights = w_mix
        )

        expect_equal(
            nobs(m_econ),
            nobs(m_lm)
        )

        expect_equal(
            df.residual(m_econ),
            df.residual(m_lm)
        )

        expect_equal(
            m_econ$rank,
            m_lm$rank
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm),
            tolerance = 1e-8
        )

        expect_equal(
            fitted(m_econ),
            fitted(m_lm),
            tolerance = 1e-8
        )

        expect_equal(
            vcov(m_econ),
            vcov(m_lm),
            tolerance = 1e-8
        )
    }
)


test_that(
    "near-collinearity produces the same rank decisions as lm",
    {

        df_aliased = mtcars %>%
            dplyr::mutate(
                wt_near = wt +
                    seq_along(wt) * 1e-14
            )

        m_econ = econ_lm(
            mpg ~ wt + wt_near + hp,
            data = df_aliased
        )

        m_lm = lm(
            mpg ~ wt + wt_near + hp,
            data = df_aliased
        )

        expect_equal(
            m_econ$rank,
            m_lm$rank
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm),
            tolerance = 1e-10
        )

        expect_equal(
            fitted(m_econ),
            fitted(m_lm),
            tolerance = 1e-10
        )

        expect_equal(
            vcov(m_econ),
            vcov(m_lm),
            tolerance = 1e-10
        )


        df_full_rank = mtcars %>%
            dplyr::mutate(
                wt_near = wt +
                    seq_along(wt) * 1e-5
            )

        m_econ = econ_lm(
            mpg ~ wt + wt_near + hp,
            data = df_full_rank
        )

        m_lm = lm(
            mpg ~ wt + wt_near + hp,
            data = df_full_rank
        )

        expect_equal(
            m_econ$rank,
            m_lm$rank
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm),
            tolerance = 1e-8
        )

        expect_equal(
            fitted(m_econ),
            fitted(m_lm),
            tolerance = 1e-8
        )

        expect_equal(
            vcov(m_econ),
            vcov(m_lm),
            tolerance = 1e-8
        )
    }
)


test_that(
    "infinite predictors and outcomes are rejected",
    {

        df_inf_x = mtcars %>%
            dplyr::mutate(
                wt = replace(
                    wt,
                    5,
                    Inf
                )
            )

        expect_error(
            econ_lm(
                mpg ~ wt + hp,
                data = df_inf_x
            ),
            "NA/NaN/Inf in 'x'"
        )


        df_inf_y = mtcars %>%
            dplyr::mutate(
                mpg = replace(
                    mpg,
                    5,
                    Inf
                )
            )

        expect_error(
            econ_lm(
                mpg ~ wt + hp,
                data = df_inf_y
            ),
            "NA/NaN/Inf in 'y'"
        )
    }
)


test_that(
    "non-finite and negative weights are rejected",
    {

        df_inf = mtcars %>%
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

        expect_error(
            econ_lm(
                mpg ~ wt + hp,
                data = df_inf,
                weights = w
            ),
            "Weights must be finite."
        )


        df_negative = mtcars %>%
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

        expect_error(
            econ_lm(
                mpg ~ wt + hp,
                data = df_negative,
                weights = w
            ),
            "Negative weights are not allowed."
        )
    }
)


test_that(
    "NaN values follow missing-data semantics",
    {

        df_x = mtcars %>%
            dplyr::mutate(
                wt = replace(
                    wt,
                    5,
                    NaN
                )
            )

        m_econ = econ_lm(
            mpg ~ wt + hp,
            data = df_x
        )

        m_lm = lm(
            mpg ~ wt + hp,
            data = df_x
        )

        expect_equal(
            nobs(m_econ),
            nobs(m_lm)
        )

        expect_equal(
            m_econ$used,
            match(
                row.names(
                    model.frame(m_lm)
                ),
                row.names(mtcars)
            )
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            vcov(m_econ),
            vcov(m_lm)
        )


        df_y = mtcars %>%
            dplyr::mutate(
                mpg = replace(
                    mpg,
                    5,
                    NaN
                )
            )

        m_econ = econ_lm(
            mpg ~ wt + hp,
            data = df_y
        )

        m_lm = lm(
            mpg ~ wt + hp,
            data = df_y
        )

        expect_equal(
            nobs(m_econ),
            nobs(m_lm)
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )
    }
)


test_that(
    "NaN weights are treated as missing weights",
    {

        df_test = mtcars %>%
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

        expect_warning(
            m_econ <- econ_lm(
                mpg ~ wt + hp,
                data = df_test,
                weights = w
            ),
            "Some observations have missingness in the weights"
        )

        m_lm = lm(
            mpg ~ wt + hp,
            data = df_test,
            weights = w
        )

        expect_equal(
            nobs(m_econ),
            nobs(m_lm)
        )

        expect_equal(
            m_econ$used,
            match(
                row.names(
                    model.frame(m_lm)
                ),
                row.names(mtcars)
            )
        )

        expect_equal(
            coef(m_econ),
            coef(m_lm)
        )

        expect_equal(
            vcov(m_econ),
            vcov(m_lm)
        )
    }
)


test_that(
    "extreme common weight scaling matches lm",
    {

        for (weight_scale in c(
            1e-300,
            1e300
        )) {

            df_test = mtcars %>%
                dplyr::mutate(
                    w = weight_scale
                )

            m_econ = suppressWarnings(
                econ_lm(
                    mpg ~ wt + hp,
                    data = df_test,
                    weights = w
                )
            )

            m_lm = lm(
                mpg ~ wt + hp,
                data = df_test,
                weights = w
            )

            expect_equal(
                nobs(m_econ),
                nobs(m_lm)
            )

            expect_equal(
                m_econ$rank,
                m_lm$rank
            )

            expect_equal(
                coef(m_econ),
                coef(m_lm),
                tolerance = 1e-10
            )

            expect_equal(
                fitted(m_econ),
                fitted(m_lm),
                tolerance = 1e-10
            )

            expect_equal(
                residuals(m_econ),
                residuals(m_lm),
                tolerance = 1e-10
            )

            expect_equal(
                vcov(m_econ),
                suppressWarnings(
                    vcov(m_lm)
                ),
                tolerance = 1e-10
            )
        }
    }
)


#### 26. API and Error Validation ####

test_that("data must be a data frame or data-frame subclass", {

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = as.matrix(mtcars)
        ),
        "'data' must be a data.frame.",
        fixed = TRUE
    )

    d = tibble::as_tibble(
        mtcars
    )

    a = econ_lm(
        mpg ~ wt + hp,
        data = d
    )

    b = lm(
        mpg ~ wt + hp,
        data = d
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )
})


test_that("invalid covariance specifications are rejected", {

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            vcov = "banana"
        )
    )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            vcov = "CR2"
        ),
        "requires a cluster variable"
    )

    d = mtcars %>%
        dplyr::mutate(
            cluster_id = rep(
                1:8,
                each = 4
            )
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = d,
            vcov = "HC1",
            cluster = cluster_id
        ),
        "is not a clustered covariance estimator"
    )
})


test_that("weights must have valid length and numeric type", {

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            weights = 1:5
        ),
        "'weights' must have one value per row of data.",
        fixed = TRUE
    )

    d_character = mtcars %>%
        dplyr::mutate(
            test_weight = rep(
                "1",
                dplyr::n()
            )
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = d_character,
            weights = test_weight
        ),
        "'weights' must be a numeric vector.",
        fixed = TRUE
    )

    d_logical = mtcars %>%
        dplyr::mutate(
            test_weight = rep(
                TRUE,
                dplyr::n()
            )
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = d_logical,
            weights = test_weight
        ),
        "'weights' must be a numeric vector.",
        fixed = TRUE
    )

    d_factor = mtcars %>%
        dplyr::mutate(
            test_weight = factor(
                rep(
                    1,
                    dplyr::n()
                )
            )
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = d_factor,
            weights = test_weight
        ),
        "'weights' must be a numeric vector.",
        fixed = TRUE
    )

    d_complex = mtcars %>%
        dplyr::mutate(
            test_weight = rep(
                1 + 1i,
                dplyr::n()
            )
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = d_complex,
            weights = test_weight
        ),
        "'weights' must be a numeric vector.",
        fixed = TRUE
    )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            weights = numeric(0)
        ),
        "'weights' must have one value per row of data.",
        fixed = TRUE
    )
})


test_that("all-zero weights are rejected", {

    d = mtcars %>%
        dplyr::mutate(
            test_weight = 0
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = d,
            weights = test_weight
        ),
        "At least one observation must have a positive weight.",
        fixed = TRUE
    )
})


test_that("empty estimation samples are rejected", {

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            subset = FALSE
        ),
        "No complete observations available for estimation.",
        fixed = TRUE
    )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            subset = logical(0)
        ),
        "No complete observations available for estimation.",
        fixed = TRUE
    )
})


test_that("missing formula variables are rejected", {

    expect_error(
        econ_lm(
            nonexistent_y ~ wt + hp,
            data = mtcars
        ),
        "object 'nonexistent_y' not found",
        fixed = TRUE
    )

    expect_error(
        econ_lm(
            mpg ~ nonexistent_x + hp,
            data = mtcars
        ),
        "object 'nonexistent_x' not found",
        fixed = TRUE
    )
})


test_that("multivariate responses are explicitly rejected", {

    expect_error(
        econ_lm(
            cbind(mpg, disp) ~ wt + hp,
            data = mtcars
        ),
        "Multivariate responses are not supported by econ_lm().",
        fixed = TRUE
    )
})


test_that("non-numeric responses are explicitly rejected", {

    d_factor = mtcars %>%
        dplyr::mutate(
            response = factor(
                mpg > median(mpg)
            )
        )

    expect_error(
        econ_lm(
            response ~ wt + hp,
            data = d_factor
        ),
        "The response variable must be numeric.",
        fixed = TRUE
    )

    d_character = mtcars %>%
        dplyr::mutate(
            response = as.character(
                mpg
            )
        )

    expect_error(
        econ_lm(
            response ~ wt + hp,
            data = d_character
        ),
        "The response variable must be numeric.",
        fixed = TRUE
    )
})


test_that("a response must be specified in the formula", {

    expect_error(
        econ_lm(
            ~ wt + hp,
            data = mtcars
        ),
        "A response variable must be specified in the formula.",
        fixed = TRUE
    )
})


test_that("numeric subsets match lm without warnings", {

    expect_no_warning(
        a <- econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            subset = c(
                1,
                3,
                5,
                10,
                20
            )
        )
    )

    b = lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = c(
            1,
            3,
            5,
            10,
            20
        )
    )

    expect_equal(
        nobs(a),
        nobs(b)
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )

    expect_equal(
        fitted(a),
        fitted(b),
        tolerance = 1e-10
    )

    expect_equal(
        vcov(a),
        vcov(b),
        tolerance = 1e-10
    )
})


test_that("logical subsets containing missing values match lm", {

    test_subset = rep(
        TRUE,
        nrow(mtcars)
    )

    test_subset[c(
        5,
        10
    )] = NA

    a = econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = test_subset
    )

    b = lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = test_subset
    )

    expect_equal(
        nobs(a),
        nobs(b)
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )

    expect_equal(
        fitted(a),
        fitted(b),
        tolerance = 1e-10
    )

    expect_equal(
        vcov(a),
        vcov(b),
        tolerance = 1e-10
    )
})


test_that("numeric subsets containing missing indices match lm", {

    test_subset = c(
        1,
        3,
        NA,
        5,
        10
    )

    expect_no_warning(
        a <- econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            subset = test_subset
        )
    )

    b = lm(
        mpg ~ wt + hp,
        data = mtcars,
        subset = test_subset
    )

    expect_equal(
        nobs(a),
        nobs(b)
    )

    expect_equal(
        coef(a),
        coef(b),
        tolerance = 1e-10
    )

    expect_equal(
        fitted(a),
        fitted(b),
        tolerance = 1e-10
    )

    expect_equal(
        vcov(a),
        vcov(b),
        tolerance = 1e-10
    )
})


test_that("invalid subset types are explicitly rejected", {

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            subset = "banana"
        ),
        "'subset' must evaluate to a logical or numeric vector.",
        fixed = TRUE
    )
})


test_that("invalid NA actions are rejected", {

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            na.action = "banana"
        )
    )
})


test_that("character and function NA actions work correctly", {

    d = mtcars %>%
        dplyr::mutate(
            hp = replace(
                hp,
                5,
                NA
            )
        )

    a_character = econ_lm(
        mpg ~ wt + hp,
        data = d,
        na.action = "na.omit"
    )

    a_function = econ_lm(
        mpg ~ wt + hp,
        data = d,
        na.action = stats::na.omit
    )

    b = lm(
        mpg ~ wt + hp,
        data = d,
        na.action = stats::na.omit
    )

    expect_equal(
        nobs(a_character),
        nobs(b)
    )

    expect_equal(
        coef(a_character),
        coef(b),
        tolerance = 1e-10
    )

    expect_equal(
        vcov(a_character),
        vcov(b),
        tolerance = 1e-10
    )

    expect_equal(
        nobs(a_function),
        nobs(b)
    )

    expect_equal(
        coef(a_function),
        coef(b),
        tolerance = 1e-10
    )

    expect_equal(
        vcov(a_function),
        vcov(b),
        tolerance = 1e-10
    )
})


test_that("cluster inputs are validated", {

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            vcov = "CR2",
            cluster = 1:5
        ),
        "'cluster' must have one value per row of data.",
        fixed = TRUE
    )

    d = mtcars %>%
        dplyr::mutate(
            cluster_a = rep(
                1:8,
                each = 4
            ),
            cluster_b = rep(
                1:4,
                each = 8
            )
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = d,
            vcov = "CR2",
            cluster = cbind(
                cluster_a,
                cluster_b
            )
        ),
        "'cluster' must have one value per row of data.",
        fixed = TRUE
    )
})


test_that("single-cluster covariance estimation is rejected", {

    d = mtcars %>%
        dplyr::mutate(
            cluster_id = 1
        )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = d,
            vcov = "CR2",
            cluster = cluster_id
        ),
        "single cluster"
    )
})


test_that("all-missing cluster IDs warn and leave no estimation sample", {

    d = mtcars %>%
        dplyr::mutate(
            cluster_id = NA_integer_
        )

    expect_warning(
        expect_error(
            econ_lm(
                mpg ~ wt + hp,
                data = d,
                vcov = "CR2",
                cluster = cluster_id
            ),
            "No complete observations available for estimation.",
            fixed = TRUE
        ),
        "Some observations have missingness in the cluster variable"
    )
})


test_that("unknown weights, subset, and cluster variables are rejected", {

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            weights = nonexistent_weights
        ),
        "object 'nonexistent_weights' not found",
        fixed = TRUE
    )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            subset = nonexistent_subset
        ),
        "object 'nonexistent_subset' not found",
        fixed = TRUE
    )

    expect_error(
        econ_lm(
            mpg ~ wt + hp,
            data = mtcars,
            vcov = "CR2",
            cluster = nonexistent_cluster
        ),
        "object 'nonexistent_cluster' not found",
        fixed = TRUE
    )
})


#### 27. ANOVA Methods ####

test_that("single-model ANOVA matches lm", {

    fit_econ = econ_lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    fit_lm = lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    expect_equal(
        anova(fit_econ),
        anova(fit_lm)
    )
})


test_that("weighted ANOVA matches lm", {

    w = seq_len(
        nrow(mtcars)
    )

    fit_econ = econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        weights = w
    )

    fit_lm = lm(
        mpg ~ wt + hp,
        data = mtcars,
        weights = w
    )

    expect_equal(
        anova(fit_econ),
        anova(fit_lm)
    )
})


test_that("ANOVA handles factors and interactions", {

    fit_econ = econ_lm(
        mpg ~ factor(cyl) * wt,
        data = mtcars
    )

    fit_lm = lm(
        mpg ~ factor(cyl) * wt,
        data = mtcars
    )

    expect_equal(
        anova(fit_econ),
        anova(fit_lm)
    )
})


test_that("ANOVA handles rank-deficient models", {

    dat = mtcars %>%
        mutate(
            wt_copy = wt
        )

    fit_econ = econ_lm(
        mpg ~ wt + wt_copy + hp,
        data = dat
    )

    fit_lm = lm(
        mpg ~ wt + wt_copy + hp,
        data = dat
    )

    expect_equal(
        anova(fit_econ),
        anova(fit_lm)
    )
})


test_that("ANOVA handles offsets", {

    fit_econ = econ_lm(
        mpg ~ wt + hp + offset(qsec),
        data = mtcars
    )

    fit_lm = lm(
        mpg ~ wt + hp + offset(qsec),
        data = mtcars
    )

    expect_equal(
        anova(fit_econ),
        anova(fit_lm)
    )
})


test_that("ANOVA handles models without intercepts", {

    fit_econ = econ_lm(
        mpg ~ 0 + wt + hp,
        data = mtcars
    )

    fit_lm = lm(
        mpg ~ 0 + wt + hp,
        data = mtcars
    )

    expect_equal(
        anova(fit_econ),
        anova(fit_lm)
    )
})


test_that("ANOVA handles zero weights", {

    w = rep(
        1,
        nrow(mtcars)
    )

    w[c(2, 7, 15)] = 0

    fit_econ = econ_lm(
        mpg ~ wt + hp,
        data = mtcars,
        weights = w
    )

    fit_lm = lm(
        mpg ~ wt + hp,
        data = mtcars,
        weights = w
    )

    expect_equal(
        anova(fit_econ),
        anova(fit_lm)
    )
})


test_that("multiple-model ANOVA matches lm", {

    econ_small = econ_lm(
        mpg ~ wt,
        data = mtcars
    )

    econ_large = econ_lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    lm_small = lm(
        mpg ~ wt,
        data = mtcars
    )

    lm_large = lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    expect_equal(
        anova(
            econ_small,
            econ_large
        ),
        anova(
            lm_small,
            lm_large
        )
    )
})


test_that("multiple-model ANOVA preserves reversed order", {

    econ_small = econ_lm(
        mpg ~ wt,
        data = mtcars
    )

    econ_large = econ_lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    lm_small = lm(
        mpg ~ wt,
        data = mtcars
    )

    lm_large = lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    expect_equal(
        anova(
            econ_large,
            econ_small
        ),
        anova(
            lm_large,
            lm_small
        )
    )
})


test_that("multiple-model ANOVA supports three models", {

    econ_0 = econ_lm(
        mpg ~ 1,
        data = mtcars
    )

    econ_1 = econ_lm(
        mpg ~ wt,
        data = mtcars
    )

    econ_2 = econ_lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    lm_0 = lm(
        mpg ~ 1,
        data = mtcars
    )

    lm_1 = lm(
        mpg ~ wt,
        data = mtcars
    )

    lm_2 = lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    expect_equal(
        anova(
            econ_0,
            econ_1,
            econ_2
        ),
        anova(
            lm_0,
            lm_1,
            lm_2
        )
    )
})


test_that("multiple-model ANOVA rejects different sample sizes", {

    dat = mtcars %>%
        mutate(
            hp = replace(
                hp,
                1,
                NA
            )
        )

    fit_1 = econ_lm(
        mpg ~ wt,
        data = dat
    )

    fit_2 = econ_lm(
        mpg ~ wt + hp,
        data = dat
    )

    expect_error(
        anova(
            fit_1,
            fit_2
        ),
        "models were not all fitted to the same size of dataset"
    )
})


test_that("multiple-model ANOVA supports scale and test arguments", {

    econ_small = econ_lm(
        mpg ~ wt,
        data = mtcars
    )

    econ_large = econ_lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    lm_small = lm(
        mpg ~ wt,
        data = mtcars
    )

    lm_large = lm(
        mpg ~ wt + hp,
        data = mtcars
    )

    expect_equal(
        anova(
            econ_small,
            econ_large,
            scale = 10
        ),
        anova(
            lm_small,
            lm_large,
            scale = 10
        )
    )

    expect_equal(
        anova(
            econ_small,
            econ_large,
            test = NULL
        ),
        anova(
            lm_small,
            lm_large,
            test = NULL
        )
    )
})