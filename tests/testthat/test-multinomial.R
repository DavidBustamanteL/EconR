###### 1. Multinomial Logit Reference Tests ####

testthat::test_that("multinomial coefficients and classical covariance match nnet", {
    set.seed(123)
    n = 600
    dat = data.frame(x1 = rnorm(n), x2 = rnorm(n), cl = rep(seq_len(60), each = 10))
    eB = 0.5 + 0.8 * dat$x1 - 0.4 * dat$x2
    eC = -0.3 + 0.3 * dat$x1 + 0.7 * dat$x2
    probs = cbind(1, exp(eB), exp(eC))
    probs = probs / rowSums(probs)
    dat$y = factor(vapply(seq_len(n), function(i) sample(c("A", "B", "C"), 1, prob = probs[i, ]), character(1)), levels = c("A", "B", "C"))

    ref = nnet::multinom(y ~ x1 + x2, data = dat, Hess = TRUE, trace = FALSE)
    fit = econ_multinomial(y ~ x1 + x2, data = dat)

    testthat::expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-5)
    testthat::expect_equal(unname(vcov(fit)), unname(vcov(ref)), tolerance = 1e-5)
    testthat::expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-5)
})


###### 2. Score-Based Covariance Tests ####

testthat::test_that("HC0 and clustered CR0 agree with sandwich", {
    set.seed(42)
    n = 900
    dat = data.frame(x1 = rnorm(n), x2 = rnorm(n), cl = rep(seq_len(90), each = 10))
    eta = cbind(0, 0.4 + dat$x1, -0.3 + dat$x2)
    P = exp(eta) / rowSums(exp(eta))
    dat$y = factor(vapply(seq_len(n), function(i) sample(c("A", "B", "C"), 1, prob = P[i, ]), character(1)), levels = c("A", "B", "C"))

    ref = nnet::multinom(y ~ x1 + x2, data = dat, Hess = TRUE, trace = FALSE)
    hc0 = econ_multinomial(y ~ x1 + x2, data = dat, vcov = "HC0")
    cr0 = econ_multinomial(y ~ x1 + x2, data = dat, vcov = "CR0", cluster = cl)

    testthat::expect_equal(
        unname(vcov(hc0)),
        unname(sandwich::sandwich(ref)),
        tolerance = 1e-4
    )

    testthat::expect_equal(
        unname(vcov(cr0)),
        unname(sandwich::vcovCL(ref, cluster = dat$cl, type = "HC0", cadjust = FALSE)),
        tolerance = 1e-4
    )

    hc1 = econ_multinomial(y ~ x1 + x2, data = dat, vcov = "HC1")
    k = length(coef(ref))

    testthat::expect_equal(
        unname(vcov(hc1)),
        unname(vcov(hc0) * n / (n - k)),
        tolerance = 1e-4
    )

    cr1s = econ_multinomial(y ~ x1 + x2, data = dat, vcov = "CR1S", cluster = cl)
    G = length(unique(dat$cl))

    testthat::expect_equal(
        unname(vcov(cr1s)),
        unname(vcov(cr0) * G * (n - 1) / ((G - 1) * (n - k))),
        tolerance = 1e-4
    )
})


###### 3. Validation Tests ####

testthat::test_that("cluster input and outcome levels are checked", {
    dat = data.frame(
        y = factor(rep(c("A", "B", "C"), each = 30)),
        x = rnorm(90)
    )

    testthat::expect_error(
        econ_multinomial(y ~ x, data = dat, vcov = "CR0"),
        "cluster"
    )

    dat$y = factor(rep(c("A", "B"), each = 45))

    testthat::expect_error(
        econ_multinomial(y ~ x, data = dat),
        "three"
    )
})



###### 4. Shared Test Data for Sample Handling ####

make_multinomial_sample = function(n = 600, seed = 2026) {

    set.seed(seed)

    dat = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        group = factor(
            rep(c("North", "South", "East"), length.out = n)
        ),
        cluster_id = rep(seq_len(n / 10), each = 10),
        w = runif(n, min = 0.5, max = 2),
        keep = rep(c(TRUE, TRUE, FALSE), length.out = n)
    )

    eta_B = 0.4 + 0.7 * dat$x1 - 0.3 * dat$x2 +
        0.4 * (dat$group == "South")

    eta_C = -0.2 + 0.2 * dat$x1 + 0.6 * dat$x2 +
        0.5 * (dat$group == "East")

    eta = cbind(0, eta_B, eta_C)

    eta = eta - apply(eta, 1, max)

    probabilities = exp(eta)
    probabilities = probabilities / rowSums(probabilities)

    dat$y = factor(
        vapply(
            seq_len(n),
            function(i) {
                sample(
                    c("A", "B", "C"),
                    size = 1,
                    prob = probabilities[i, ]
                )
            },
            character(1)
        ),
        levels = c("A", "B", "C")
    )

    dat
}


###### 5. Weighted Multinomial Estimation ####

testthat::test_that("weighted multinomial matches nnet", {

    dat = make_multinomial_sample()

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = dat,
        weights = w,
        Hess = TRUE,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        weights = w
    )

    testthat::expect_equal(
        unname(coef(fit)),
        unname(coef(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        unname(vcov(fit)),
        unname(vcov(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        as.numeric(logLik(fit)),
        as.numeric(logLik(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        nobs(fit),
        nrow(dat)
    )
})


###### 6. Logical Subset Estimation ####

testthat::test_that("logical subsets match nnet", {

    dat = make_multinomial_sample()

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = dat,
        subset = keep,
        Hess = TRUE,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        subset = keep
    )

    testthat::expect_equal(
        unname(coef(fit)),
        unname(coef(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        unname(vcov(fit)),
        unname(vcov(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        nobs(fit),
        sum(dat$keep)
    )

    testthat::expect_equal(
        fit$used,
        which(dat$keep)
    )
})


###### 7. Weighted Subset Estimation ####

testthat::test_that("weights and subsets work together", {

    dat = make_multinomial_sample()

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
        subset = keep
    )

    testthat::expect_equal(
        unname(coef(fit)),
        unname(coef(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        unname(vcov(fit)),
        unname(vcov(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        as.numeric(logLik(fit)),
        as.numeric(logLik(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        fit$econ_weights,
        dat$w[dat$keep]
    )
})


###### 8. Missing Values and Sample Alignment ####

testthat::test_that("missing values match nnet sample selection", {

    dat = make_multinomial_sample()

    dat$x1[c(5, 20, 100)] = NA
    dat$x2[c(10, 30)] = NA

    complete = stats::complete.cases(
        dat[, c("y", "x1", "x2")]
    )

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = dat,
        na.action = stats::na.omit,
        Hess = TRUE,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        na.action = stats::na.omit
    )

    testthat::expect_equal(
        unname(coef(fit)),
        unname(coef(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        unname(vcov(fit)),
        unname(vcov(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        nobs(fit),
        sum(complete)
    )

    testthat::expect_equal(
        fit$used,
        which(complete)
    )

    testthat::expect_equal(
        fit$omitted,
        which(!complete)
    )
})


###### 9. Cluster Alignment After Missing Values ####

testthat::test_that("cluster scores use the estimation sample", {

    dat = make_multinomial_sample()

    dat$x1[c(5, 20, 100)] = NA
    dat$x2[c(10, 30)] = NA

    complete = stats::complete.cases(
        dat[, c("y", "x1", "x2")]
    )

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = dat,
        na.action = stats::na.omit,
        Hess = TRUE,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        vcov = "CR0",
        cluster = cluster_id
    )

    V_ref = sandwich::vcovCL(
        ref,
        cluster = dat$cluster_id[complete],
        type = "HC0",
        cadjust = FALSE
    )

    testthat::expect_equal(
        unname(vcov(fit)),
        unname(V_ref),
        tolerance = 1e-4
    )

    testthat::expect_equal(
        fit$cluster,
        dat$cluster_id[complete]
    )
})


###### 10. Factor Predictor Estimation ####

testthat::test_that("factor predictors match nnet", {

    dat = make_multinomial_sample()

    ref = nnet::multinom(
        y ~ x1 + group,
        data = dat,
        Hess = TRUE,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + group,
        data = dat
    )

    testthat::expect_equal(
        unname(coef(fit)),
        unname(coef(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        unname(vcov(fit)),
        unname(vcov(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        colnames(coef(fit)),
        colnames(coef(ref))
    )
})


###### 11. Invalid and Missing Weight Handling ####

testthat::test_that("invalid weights are rejected and missing weights are dropped", {

    dat = make_multinomial_sample()

    ###### 1. Reject negative weights ####

    dat$w[1] = -1

    testthat::expect_error(
        econ_multinomial(
            y ~ x1 + x2,
            data = dat,
            weights = w
        ),
        "weights"
    )

    ###### 2. Reject infinite weights ####

    dat$w[1] = Inf

    testthat::expect_error(
        econ_multinomial(
            y ~ x1 + x2,
            data = dat,
            weights = w
        ),
        "weights"
    )

    ###### 3. Exclude missing weights ####

    dat$w[1] = NA_real_

    dat$w[20] = NA_real_

    complete = !is.na(dat$w)

    testthat::expect_warning(
    {
        fit = econ_multinomial(
            y ~ x1 + x2,
            data = dat,
            weights = w
        )
    },
    "missingness in the weights"
    )

    ###### 4. Validate estimation sample ####

    testthat::expect_equal(
        nobs(fit),
        sum(complete)
    )

    testthat::expect_equal(
        fit$used,
        which(complete)
    )

    testthat::expect_equal(
        fit$omitted,
        which(!complete)
    )

    testthat::expect_equal(
        unname(fit$econ_weights),
        unname(dat$w[complete])
    )

    ###### 5. Compare with nnet ####

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = dat,
        weights = w,
        Hess = TRUE,
        trace = FALSE
    )

    testthat::expect_equal(
        unname(coef(fit)),
        unname(coef(ref)),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        unname(vcov(fit)),
        unname(vcov(ref)),
        tolerance = 1e-5
    )

})



###### 12. Multinomial Link Selection ####

testthat::test_that("multinomial logit is the default link", {

    dat = make_multinomial_sample()

    model_default = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    model_logit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        link = "logit"
    )

    testthat::expect_identical(
        model_default$link,
        "logit"
    )

    testthat::expect_identical(
        model_logit$link,
        "logit"
    )

    testthat::expect_equal(
        coef(model_default),
        coef(model_logit),
        tolerance = 1e-10
    )

    testthat::expect_equal(
        vcov(model_default),
        vcov(model_logit),
        tolerance = 1e-10
    )
})


###### 13. Unsupported Multinomial Probit ####

testthat::test_that("multinomial link selection rejects invalid links", {

    dat = make_multinomial_sample(n = 90, seed = 1301)

    testthat::expect_error(
        econ_multinomial(
            y ~ x1 + x2,
            data = dat,
            link = "invalid"
        ),
        "arg"
    )
})



###### 14. EconR Summary Uses Selected Covariance ####

testthat::test_that("summary uses the selected covariance estimator", {

    dat = make_multinomial_sample()

    model = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        link = "logit",
        vcov = "HC1"
    )

    s = summary(model)

    testthat::expect_s3_class(
        s,
        "summary.econ_multinomial"
    )

    testthat::expect_identical(
        s$vcov_type,
        "HC1"
    )

    expected_se = matrix(
        sqrt(diag(vcov(model))),
        nrow = nrow(coef(model)),
        byrow = TRUE,
        dimnames = dimnames(coef(model))
    )

    testthat::expect_equal(
        s$standard.errors,
        expected_se,
        tolerance = 1e-10
    )
})



###### 15. Cluster-Robust Covariance Corrections ####

testthat::test_that("all cluster corrections match independent references", {

    dat = make_multinomial_sample(
        n = 1000,
        seed = 1234
    )

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = dat,
        Hess = TRUE,
        trace = FALSE
    )

    n = nrow(dat)
    k = length(coef(ref))
    G = length(unique(dat$cluster_id))

    V_CR0 = sandwich::vcovCL(
        ref,
        cluster = dat$cluster_id,
        type = "HC0",
        cadjust = FALSE
    )

    factors = c(
        CR0 = 1,
        CR1 = (n - 1) / (n - k),
        CR1p = G / (G - k),
        CR1S = G * (n - 1) / ((G - 1) * (n - k)),
        stata = G * (n - 1) / ((G - 1) * (n - k))
    )

    for (type in names(factors)) {

        fit = econ_multinomial(
            y ~ x1 + x2,
            data = dat,
            vcov = type,
            cluster = cluster_id
        )

        testthat::expect_equal(
            unname(vcov(fit)),
            unname(V_CR0 * factors[[type]]),
            tolerance = 1e-4,
            info = paste("Covariance mismatch:", type)
        )
    }
})


###### 16. Clustered Inference and Degrees of Freedom ####

testthat::test_that("clustered inference uses G minus 1 degrees of freedom", {

    dat = make_multinomial_sample(
        n = 600,
        seed = 1235
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        vcov = "stata",
        cluster = cluster_id
    )

    s = summary(fit)

    G = length(unique(dat$cluster_id))

    testthat::expect_equal(
        s$df,
        G - 1
    )

    expected_se = matrix(
        sqrt(diag(vcov(fit))),
        nrow = nrow(coef(fit)),
        byrow = TRUE,
        dimnames = dimnames(coef(fit))
    )

    testthat::expect_equal(
        s$standard.errors,
        expected_se,
        tolerance = 1e-10
    )

    expected_p = 2 * stats::pt(
        -abs(s$statistic),
        df = G - 1
    )

    testthat::expect_equal(
        s$p.value,
        expected_p,
        tolerance = 1e-10
    )
})


###### 17. Nonclustered Inference ####

testthat::test_that("classical and HC inference use normal approximation", {

    dat = make_multinomial_sample(
        n = 600,
        seed = 1236
    )

    for (type in c("classical", "HC0", "HC1")) {

        fit = econ_multinomial(
            y ~ x1 + x2,
            data = dat,
            vcov = type
        )

        s = summary(fit)

        testthat::expect_true(
            is.infinite(s$df)
        )

        expected_p = 2 * stats::pnorm(
            -abs(s$statistic)
        )

        testthat::expect_equal(
            s$p.value,
            expected_p,
            tolerance = 1e-10
        )
    }
})



###### 18. Weighted and Subset Clustered Covariance ####

testthat::test_that("clustered covariance handles weights and subsets", {

    dat = make_multinomial_sample(
        n = 1000,
        seed = 1238
    )

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

    ###### Independent Weighted Score Calculation ####

    selected = which(dat$keep)

    X = stats::model.matrix(
        ~ x1 + x2,
        data = dat[selected, , drop = FALSE]
    )

    P = predict(ref, type = "probs")

    y = dat$y[selected]
    w = dat$w[selected]

    S_expected = do.call(
        cbind,
        lapply(levels(y)[-1], function(category) {

            residual = as.numeric(y == category) -
                P[, category]

            X * (w * residual)
        })
    )

    ###### Independent Clustered Sandwich Calculation ####

    cluster = dat$cluster_id[selected]

    U = rowsum(
        S_expected,
        group = cluster,
        reorder = TRUE
    )

    V_classical = stats::vcov(ref)

    V_expected = V_classical %*%
        crossprod(U) %*%
        V_classical

    ###### Validate Observation Scores ####

    testthat::expect_equal(
        unname(sandwich::estfun(fit)),
        unname(S_expected),
        tolerance = 1e-7
    )

    ###### Validate Cluster-Robust Covariance ####

    testthat::expect_equal(
        unname(vcov(fit)),
        unname(V_expected),
        tolerance = 1e-7
    )

    ###### Validate Estimation Sample and Clusters ####

    testthat::expect_equal(
        nobs(fit),
        length(selected)
    )

    testthat::expect_equal(
        fit$cluster,
        cluster
    )

    testthat::expect_equal(
        as.numeric(fit$econ_weights),
        w
    )

    ###### Validate Weighted Score Equations ####

    testthat::expect_equal(
        as.numeric(colSums(sandwich::estfun(fit))),
        rep(0, ncol(S_expected)),
        tolerance = 0.02
    )
})



###### 19. In-Sample Predicted Probabilities ####

testthat::test_that("multinomial probabilities match nnet", {

    dat = make_multinomial_sample(n = 1000, seed = 2101)

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = dat,
        Hess = TRUE,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    P_ref = predict(ref, type = "probs")
    P_fit = predict(fit, type = "probs")

    testthat::expect_equal(
        unname(P_fit),
        unname(P_ref),
        tolerance = 1e-8
    )

    testthat::expect_equal(
        unname(rowSums(P_fit)),
        rep(1, nrow(dat)),
        tolerance = 1e-12
    )

    testthat::expect_true(all(is.finite(P_fit)))
    testthat::expect_true(all(P_fit >= 0))
    testthat::expect_true(all(P_fit <= 1))

    testthat::expect_identical(
        colnames(P_fit),
        levels(dat$y)
    )
})


###### 20. Predicted Outcome Categories ####

testthat::test_that("predicted classes match probability maxima", {

    dat = make_multinomial_sample(n = 1000, seed = 2102)

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = dat,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    class_ref = predict(ref, type = "class")
    class_fit = predict(fit, type = "class")

    P = predict(fit, type = "probs")

    expected = colnames(P)[
        max.col(P, ties.method = "first")
    ]

    testthat::expect_equal(
        as.character(class_fit),
        as.character(class_ref)
    )

    testthat::expect_equal(
        as.character(class_fit),
        expected
    )

    testthat::expect_identical(
        levels(class_fit),
        levels(dat$y)
    )
})


###### 21. Out-of-Sample Prediction ####

testthat::test_that("newdata predictions match nnet", {

    dat = make_multinomial_sample(n = 1000, seed = 2103)

    training = dat[1:800, ]
    testing = dat[801:1000, ]

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = training,
        Hess = TRUE,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = training
    )

    P_ref = predict(
        ref,
        newdata = testing,
        type = "probs"
    )

    P_fit = predict(
        fit,
        newdata = testing,
        type = "probs"
    )

    testthat::expect_equal(
        unname(P_fit),
        unname(P_ref),
        tolerance = 1e-8
    )

    testthat::expect_equal(
        dim(P_fit),
        c(200L, 3L)
    )

    testthat::expect_equal(
        unname(rowSums(P_fit)),
        rep(1, 200),
        tolerance = 1e-12
    )

    testthat::expect_equal(
        as.character(predict(fit, testing, type = "class")),
        as.character(predict(ref, testing, type = "class"))
    )
})


###### 22. Factor Predictors and New Levels ####

testthat::test_that("factor predictions match nnet", {

    dat = make_multinomial_sample(n = 1000, seed = 2104)

    training = dat[1:800, ]
    testing = dat[801:1000, ]

    ref = nnet::multinom(
        y ~ x1 + group,
        data = training,
        Hess = TRUE,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + group,
        data = training
    )

    testthat::expect_equal(
        unname(predict(fit, testing, type = "probs")),
        unname(predict(ref, testing, type = "probs")),
        tolerance = 1e-8
    )

    invalid = testing[1:5, ]

    invalid$group = as.character(invalid$group)
    invalid$group[1] = "Unknown"
    invalid$group = factor(invalid$group)

    testthat::expect_error(
        predict(fit, newdata = invalid, type = "probs")
    )
})


###### 23. Reference Category Invariance ####

testthat::test_that("reference category changes coefficients but not probabilities", {

    dat = make_multinomial_sample(n = 1000, seed = 2105)

    fit_A = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    dat$y = stats::relevel(dat$y, ref = "C")

    fit_C = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    P_A = predict(fit_A, type = "probs")
    P_C = predict(fit_C, type = "probs")

    P_C = P_C[, colnames(P_A), drop = FALSE]

    testthat::expect_equal(
        unname(P_A),
        unname(P_C),
        tolerance = 1e-4
    )

    testthat::expect_equal(
        as.numeric(logLik(fit_A)),
        as.numeric(logLik(fit_C)),
        tolerance = 1e-4
    )

    testthat::expect_identical(
        rownames(coef(fit_C)),
        c("A", "B")
    )

    testthat::expect_identical(
        colnames(predict(fit_C, type = "probs")),
        c("C", "A", "B")
    )
})


###### 24. Missing Values in Prediction Data ####

testthat::test_that("missing newdata predictions match nnet", {

    dat = make_multinomial_sample(n = 1000, seed = 2106)

    training = dat[1:800, ]
    testing = dat[801:1000, ]

    testing$x1[c(2, 5)] = NA_real_
    testing$x2[10] = NA_real_

    ref = nnet::multinom(
        y ~ x1 + x2,
        data = training,
        trace = FALSE
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = training
    )

    P_ref = predict(
        ref,
        newdata = testing,
        type = "probs"
    )

    P_fit = predict(
        fit,
        newdata = testing,
        type = "probs"
    )

    testthat::expect_equal(
        P_fit,
        P_ref,
        tolerance = 1e-8
    )

    missing_rows = which(
        !stats::complete.cases(testing[, c("x1", "x2")])
    )

    testthat::expect_true(
        all(is.na(P_fit[missing_rows, , drop = FALSE]))
    )

    testthat::expect_equal(
        unname(rowSums(P_fit[-missing_rows, , drop = FALSE])),
        rep(1, nrow(testing) - length(missing_rows)),
        tolerance = 1e-12
    )
})


###### 25. Scores Under Alternative Reference Categories ####

testthat::test_that("scores match analytical formula after changing baseline", {

    dat = make_multinomial_sample(n = 1000, seed = 2107)

    dat$y = stats::relevel(dat$y, ref = "C")

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    X = stats::model.matrix(
        ~ x1 + x2,
        data = dat
    )

    P = predict(fit, type = "probs")

    nonbase = levels(dat$y)[-1]

    S_expected = do.call(
        cbind,
        lapply(nonbase, function(category) {

            residual = as.numeric(dat$y == category) -
                P[, category]

            X * residual
        })
    )

    S_fit = sandwich::estfun(fit)

    testthat::expect_equal(
        unname(S_fit),
        unname(S_expected),
        tolerance = 1e-8
    )

    testthat::expect_equal(
        ncol(S_fit),
        length(coef(fit))
    )

    testthat::expect_identical(
        colnames(S_fit),
        colnames(vcov(fit))
    )
})



###### 26. Classical Wald Confidence Intervals ####

testthat::test_that("classical multinomial confidence intervals are correct", {

    dat = make_multinomial_sample(n = 1000, seed = 2201)

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        vcov = "classical"
    )

    V = vcov(fit)
    estimates = as.vector(t(coef(fit)))
    names(estimates) = colnames(V)

    critical = stats::qnorm(0.975)

    expected = cbind(
        estimates - critical * sqrt(diag(V)),
        estimates + critical * sqrt(diag(V))
    )

    actual = confint(fit)

    testthat::expect_equal(
        unname(actual),
        unname(expected),
        tolerance = 1e-10
    )

    testthat::expect_equal(dim(actual), c(6L, 2L))
    testthat::expect_identical(rownames(actual), colnames(V))
})


###### 27. Clustered Wald Confidence Intervals ####

testthat::test_that("clustered intervals use G minus 1 degrees of freedom", {

    dat = make_multinomial_sample(n = 1000, seed = 2202)

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        vcov = "CR1S",
        cluster = cluster_id
    )

    V = vcov(fit)
    estimates = as.vector(t(coef(fit)))
    names(estimates) = colnames(V)

    G = length(unique(fit$cluster))
    critical = stats::qt(0.975, df = G - 1L)

    expected = cbind(
        estimates - critical * sqrt(diag(V)),
        estimates + critical * sqrt(diag(V))
    )

    testthat::expect_equal(
        unname(confint(fit)),
        unname(expected),
        tolerance = 1e-10
    )
})


###### 28. Confidence Interval Arguments ####

testthat::test_that("confidence interval arguments are validated", {

    dat = make_multinomial_sample(n = 600, seed = 2203)

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    full = confint(fit)
    selected = rownames(full)[c(2, 5)]

    testthat::expect_equal(
        confint(fit, parm = selected),
        full[selected, , drop = FALSE]
    )

    testthat::expect_equal(
        nrow(confint(fit, level = 0.90)),
        nrow(full)
    )

    testthat::expect_error(confint(fit, level = 0))
    testthat::expect_error(confint(fit, level = 1))
    testthat::expect_error(confint(fit, level = NA_real_))
    testthat::expect_error(confint(fit, parm = "nonexistent"))
    testthat::expect_error(confint(fit, parm = c(selected[1], selected[1])))
})



###### 29. Classical Joint Wald Test ####

testthat::test_that("classical multinomial joint Wald test is correct", {

    dat = make_multinomial_sample(n = 1000, seed = 2301)

    fit = econ_multinomial(y ~ x1 + x2, data = dat)

    s = summary(fit)

    V = vcov(fit)
    beta = as.vector(t(coef(fit)))
    names(beta) = colnames(V)

    selected = names(beta)[
    !endsWith(names(beta), ":(Intercept)")
    ]

    b = beta[selected]
    V_test = V[selected, selected, drop = FALSE]

    W = as.numeric(crossprod(b, solve(V_test, b)))
    q = length(b)

    testthat::expect_equal(s$wald.statistic, W)
    testthat::expect_equal(s$wald.df1, q)
    testthat::expect_identical(s$wald.distribution, "Chisq")

    testthat::expect_equal(
        s$wald.p.value,
        pchisq(W, df = q, lower.tail = FALSE)
    )
})


###### 30. Clustered Joint Wald Test ####

testthat::test_that("clustered multinomial Wald test uses F distribution", {

    dat = make_multinomial_sample(n = 1000, seed = 2302)

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        vcov = "stata",
        cluster = cluster_id
    )

    s = summary(fit)

    V = vcov(fit)
    beta = as.vector(t(coef(fit)))
    names(beta) = colnames(V)

    selected = names(beta)[
    !endsWith(names(beta), ":(Intercept)")
    ]

    b = beta[selected]
    V_test = V[selected, selected, drop = FALSE]

    W = as.numeric(crossprod(b, solve(V_test, b)))
    q = length(b)
    G = length(unique(fit$cluster))

    testthat::expect_identical(s$wald.distribution, "F")
    testthat::expect_equal(s$wald.statistic, W / q)
    testthat::expect_equal(s$wald.df1, q)
    testthat::expect_equal(s$wald.df2, G - 1L)

    testthat::expect_equal(
        s$wald.p.value,
        pf(W / q, q, G - 1L, lower.tail = FALSE)
    )
})


###### 31. Joint Wald Tests Across Covariance Types ####

testthat::test_that("joint Wald tests use the selected covariance", {

    dat = make_multinomial_sample(n = 1000, seed = 2303)

    for (type in c(
        "classical", "HC0", "HC1",
        "CR0", "CR1", "CR1p", "CR1S", "stata"
    )) {

        fit = econ_multinomial(
            y ~ x1 + x2,
            data = dat,
            vcov = type,
            cluster = cluster_id
        )

        s = summary(fit)

        testthat::expect_true(
            is.finite(s$wald.statistic)
        )

        testthat::expect_true(
            is.finite(s$wald.p.value)
        )

        testthat::expect_true(
            s$wald.p.value >= 0 &&
                s$wald.p.value <= 1
        )

        testthat::expect_identical(
            s$wald.distribution,
            if (type %in% c("classical", "HC0", "HC1")) {
                "Chisq"
            } else {
                "F"
            }
        )
    }
})



###### 32. Multinomial Average Marginal Effects ####

testthat::test_that("multinomial AMEs satisfy adding-up", {

    dat = make_multinomial_sample(n = 1000, seed = 2401)

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    ame = marginal_effects.econ_multinomial(fit)

    testthat::expect_equal(nrow(ame), 6L)

    for (v in c("x1", "x2")) {

        testthat::expect_equal(
            sum(ame$AME[ame$variable == v]),
            0,
            tolerance = 1e-10
        )
    }
})


###### 33. AMEs Against Numerical Probability Derivatives ####

testthat::test_that("multinomial AMEs match numerical derivatives", {

    dat = make_multinomial_sample(n = 1000, seed = 2402)

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    ame = marginal_effects.econ_multinomial(
        fit,
        variables = "x1"
    )

    h = 1e-5

    plus = dat
    minus = dat

    plus$x1 = plus$x1 + h
    minus$x1 = minus$x1 - h

    P_plus = predict(fit, newdata = plus, type = "probs")
    P_minus = predict(fit, newdata = minus, type = "probs")

    numerical = colMeans(
        (P_plus - P_minus) / (2 * h)
    )

    testthat::expect_equal(
        unname(ame$AME),
        unname(numerical[ame$outcome]),
        tolerance = 1e-6
    )
})


###### 34. Clustered AME Inference ####

testthat::test_that("clustered AMEs use cluster-adjusted inference", {

    dat = make_multinomial_sample(n = 1000, seed = 2403)

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        vcov = "stata",
        cluster = cluster_id
    )

    ame = marginal_effects.econ_multinomial(
        fit,
        variables = "x1"
    )

    G = length(unique(fit$cluster))

    testthat::expect_equal(
        ame$p.value,
        2 * pt(
            -abs(ame$statistic),
            df = G - 1L
        )
    )

    testthat::expect_true(
        all(is.finite(ame$std.error))
    )

    testthat::expect_true(
        all(ame$std.error >= 0)
    )
})


###### 35. Weighted Average Marginal Effects ####

testthat::test_that("weighted multinomial AMEs are computed correctly", {

    dat = make_multinomial_sample(n = 1000, seed = 2405)

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat,
        weights = w
    )

    ame_unweighted = marginal_effects.econ_multinomial(
        fit,
        variables = "x1",
        weighted = FALSE
    )

    ame_weighted = marginal_effects.econ_multinomial(
        fit,
        variables = "x1",
        weighted = TRUE
    )

    testthat::expect_equal(
        nrow(ame_weighted),
        3L
    )

    testthat::expect_equal(
        sum(ame_weighted$AME),
        0,
        tolerance = 1e-10
    )

    testthat::expect_true(
        all(is.finite(ame_weighted$std.error))
    )

    testthat::expect_true(
        any(abs(
            ame_weighted$AME - ame_unweighted$AME
        ) > 1e-8)
    )

    ###### Independent weighted derivative check ####

    h = 1e-5

    plus = dat
    minus = dat

    plus$x1 = plus$x1 + h
    minus$x1 = minus$x1 - h

    P_plus = predict(fit, newdata = plus, type = "probs")
    P_minus = predict(fit, newdata = minus, type = "probs")

    derivatives = (P_plus - P_minus) / (2 * h)

    reference = colSums(
        derivatives * dat$w
    ) / sum(dat$w)

    testthat::expect_equal(
        unname(ame_weighted$AME),
        unname(reference[ame_weighted$outcome]),
        tolerance = 1e-6
    )
})



###### 36. Categorical Average Marginal Effects ####

testthat::test_that("categorical AMEs match discrete probability changes", {

    dat = make_multinomial_sample(n = 1000, seed = 2501)

    fit = econ_multinomial(
        y ~ x1 + group,
        data = dat
    )

    ame = marginal_effects.econ_multinomial(
        fit,
        variables = "group"
    )

    reference = levels(dat$group)[1L]
    targets = levels(dat$group)[-1L]

    testthat::expect_equal(
        nrow(ame),
        length(targets) * 3L
    )

    for (target in targets) {

        high = dat
        low = dat

        high$group = factor(
            rep(target, nrow(dat)),
            levels = levels(dat$group)
        )

        low$group = factor(
            rep(reference, nrow(dat)),
            levels = levels(dat$group)
        )

        P_high = predict(fit, newdata = high, type = "probs")
        P_low = predict(fit, newdata = low, type = "probs")

        expected = colMeans(P_high - P_low)

        contrast_name = paste0(target, " - ", reference)

        actual = ame[
            ame$contrast == contrast_name,
            ,
            drop = FALSE
        ]

        testthat::expect_equal(
            unname(actual$AME),
            unname(expected[actual$outcome]),
            tolerance = 1e-8
        )

        testthat::expect_equal(
            sum(actual$AME),
            0,
            tolerance = 1e-10
        )
    }

    testthat::expect_true(
        all(is.finite(ame$std.error))
    )
})


testthat::test_that("weighted categorical AMEs match weighted contrasts", {

    dat = make_multinomial_sample(n = 1000, seed = 2502)

    fit = econ_multinomial(
        y ~ x1 + group,
        data = dat,
        weights = w
    )

    ame = marginal_effects.econ_multinomial(
        fit,
        variables = "group",
        weighted = TRUE
    )

    reference = levels(dat$group)[1L]

    for (target in levels(dat$group)[-1L]) {

        high = dat
        low = dat

        high$group = factor(
            rep(target, nrow(dat)),
            levels = levels(dat$group)
        )

        low$group = factor(
            rep(reference, nrow(dat)),
            levels = levels(dat$group)
        )

        difference =
            predict(fit, newdata = high, type = "probs") -
            predict(fit, newdata = low, type = "probs")

        expected = colSums(
            difference * dat$w
        ) / sum(dat$w)

        actual = ame[
            ame$contrast == paste0(target, " - ", reference),
            ,
            drop = FALSE
        ]

        testthat::expect_equal(
            unname(actual$AME),
            unname(expected[actual$outcome]),
            tolerance = 1e-8
        )
    }
})


testthat::test_that("continuous and categorical AMEs work together", {

    dat = make_multinomial_sample(n = 1000, seed = 2503)

    fit = econ_multinomial(
        y ~ x1 + group,
        data = dat
    )

    ame = marginal_effects.econ_multinomial(fit)

    testthat::expect_equal(
        unique(ame$variable),
        c("x1", "group")
    )

    testthat::expect_equal(
        nrow(ame),
        9L
    )

    testthat::expect_equal(
        sum(ame$AME[ame$variable == "x1"]),
        0,
        tolerance = 1e-10
    )
})



###### 37. Multinomial Augment Tests ####

testthat::test_that("multinomial augment matches fitted predictions", {

    set.seed(2406)

    n = 300

    dat = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        region = factor(
            sample(c("East", "North", "South"), n, replace = TRUE)
        )
    )

    eta_B = 0.4 + 0.7 * dat$x1 - 0.3 * dat$x2
    eta_C = -0.2 + 0.2 * dat$x1 + 0.6 * dat$x2

    eta = cbind(A = 0, B = eta_B, C = eta_C)
    probabilities = exp(eta) / rowSums(exp(eta))

    dat$y = factor(
        vapply(seq_len(n), function(i) {
            sample(
                c("A", "B", "C"),
                size = 1,
                prob = probabilities[i, ]
            )
        }, character(1)),
        levels = c("A", "B", "C")
    )

    fit = econ_multinomial(
        y ~ x1 + x2 + region,
        data = dat
    )

    aug = broom::augment(fit)
    predicted_probs = predict(fit, type = "probs")
    predicted_class = predict(fit, type = "class")

    testthat::expect_s3_class(aug, "tbl_df")
    testthat::expect_equal(nrow(aug), nobs(fit))

    testthat::expect_identical(
        as.character(aug$.fitted),
        as.character(predicted_class)
    )

    testthat::expect_equal(
        unname(as.matrix(
            aug[, paste0(".fitted_", colnames(predicted_probs))]
        )),
        unname(as.matrix(predicted_probs)),
        tolerance = 1e-10
    )

    testthat::expect_equal(
        unname(rowSums(predicted_probs)),
        rep(1, nobs(fit)),
        tolerance = 1e-10
    )

    testthat::expect_false("(weights)" %in% names(aug))
})

testthat::test_that("multinomial augment supports newdata", {

    set.seed(2407)

    n = 300

    dat = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        y = factor(
            sample(c("A", "B", "C"), n, replace = TRUE),
            levels = c("A", "B", "C")
        )
    )

    fit = econ_multinomial(y ~ x1 + x2, data = dat)

    new_dat = dat[1:10, c("x1", "x2")]

    aug_new = broom::augment(fit, newdata = new_dat)

    predicted_probs = predict(
        fit,
        newdata = new_dat,
        type = "probs"
    )

    predicted_class = predict(
        fit,
        newdata = new_dat,
        type = "class"
    )

    testthat::expect_equal(nrow(aug_new), 10L)

    testthat::expect_identical(
        as.character(aug_new$.fitted),
        as.character(predicted_class)
    )

    testthat::expect_equal(
        unname(as.matrix(
            aug_new[, paste0(".fitted_", colnames(predicted_probs))]
        )),
        unname(as.matrix(predicted_probs)),
        tolerance = 1e-10
    )

    testthat::expect_error(
        broom::augment(
            fit,
            data = dat,
            newdata = new_dat
        ),
        "either"
    )
})



###### 38. Multinomial McFadden Pseudo-R-Squared Tests ####

testthat::test_that("multinomial pseudo-R-squared matches reference calculation", {

    set.seed(2409)

    n = 400

    dat = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        y = factor(
            sample(c("A", "B", "C"), n, replace = TRUE),
            levels = c("A", "B", "C")
        )
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    result = broom::glance(fit)

    counts = as.numeric(table(dat$y))
    probabilities = counts / sum(counts)

    null_loglik = sum(counts * log(probabilities))

    reference_r2 = 1 -
        as.numeric(stats::logLik(fit)) / null_loglik

    testthat::expect_equal(
        result$pseudo.r.squared,
        reference_r2,
        tolerance = 1e-12
    )

    testthat::expect_true(
        all(c(
            "edf", "deviance", "AIC", "nobs",
            "pseudo.r.squared"
        ) %in% names(result))
    )

    testthat::expect_equal(
        result$nobs,
        stats::nobs(fit)
    )

    testthat::expect_equal(
        result$AIC,
        stats::AIC(fit)
    )
})


testthat::test_that("multinomial pseudo-R-squared respects estimation weights", {

    set.seed(2410)

    n = 400

    dat = data.frame(
        x1 = rnorm(n),
        y = factor(
            sample(c("A", "B", "C"), n, replace = TRUE),
            levels = c("A", "B", "C")
        ),
        w = sample(1:4, n, replace = TRUE)
    )

    fit = econ_multinomial(
        y ~ x1,
        data = dat,
        weights = w
    )

    result = broom::glance(fit)

    counts = as.numeric(
        tapply(dat$w, dat$y, sum)
    )

    probabilities = counts / sum(counts)

    null_loglik = sum(
        counts * log(probabilities)
    )

    reference_r2 = 1 -
        as.numeric(stats::logLik(fit)) / null_loglik

    testthat::expect_equal(
        result$pseudo.r.squared,
        reference_r2,
        tolerance = 1e-12
    )

    testthat::expect_true(
        is.finite(result$pseudo.r.squared)
    )
})



###### 39. Multinomial Summary Tests ####

testthat::test_that("multinomial summary reports correct model information", {

    set.seed(2411)

    n = 500

    dat = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        y = factor(
            sample(c("A", "B", "C"), n, replace = TRUE),
            levels = c("A", "B", "C")
        )
    )

    fit = econ_multinomial(
        y ~ x1 + x2,
        data = dat
    )

    s = summary(fit)

    testthat::expect_identical(
        s$call,
        fit$call
    )

    testthat::expect_identical(
        fit$call[[1]],
        as.name("econ_multinomial")
    )

    testthat::expect_equal(
        s$pseudo.r.squared,
        broom::glance(fit)$pseudo.r.squared,
        tolerance = 1e-12
    )

    testthat::expect_equal(
        s$standard.errors^2,
        matrix(
            diag(vcov(fit)),
            nrow = nrow(s$standard.errors),
            byrow = TRUE,
            dimnames = dimnames(s$standard.errors)
        ),
        tolerance = 1e-10
    )

    testthat::expect_identical(
        s$extra.significance,
        FALSE
    )

    output = capture.output(print(s))

    testthat::expect_true(
        any(grepl("McFadden pseudo-R-squared:", output, fixed = TRUE))
    )

    testthat::expect_true(
        any(grepl("econ_multinomial(", output, fixed = TRUE))
    )
})


testthat::test_that("multinomial optional significance does not change inference", {

    set.seed(2412)

    n = 500

    dat = data.frame(
        x1 = rnorm(n),
        y = factor(
            sample(c("A", "B", "C"), n, replace = TRUE),
            levels = c("A", "B", "C")
        )
    )

    fit = econ_multinomial(y ~ x1, data = dat)

    s = summary(fit)
    s_extra = summary(fit, extra.significance = TRUE)

    testthat::expect_identical(
        s_extra$extra.significance,
        TRUE
    )

    testthat::expect_identical(
        s$coefficients,
        s_extra$coefficients
    )

    testthat::expect_identical(
        s$p.value,
        s_extra$p.value
    )

    testthat::expect_identical(
        s$wald.statistic,
        s_extra$wald.statistic
    )

    output = capture.output(print(s_extra))

    testthat::expect_true(
        any(grepl("+ p < 0.001", output, fixed = TRUE))
    )

    testthat::expect_error(
        summary(fit, extra.significance = NA),
        "TRUE or FALSE"
    )
})



###### 40. Multinomial Tidy Tests ####

testthat::test_that("multinomial tidy matches classical and clustered inference", {

    set.seed(2414)

    n = 500

    dat = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        cluster_id = sample(seq_len(40), n, replace = TRUE),
        y = factor(
            sample(c("A", "B", "C"), n, replace = TRUE),
            levels = c("A", "B", "C")
        )
    )

    for (type in c("classical", "HC1", "stata")) {

        fit = if (type == "stata") {
            econ_multinomial(
                y ~ x1 + x2,
                data = dat,
                vcov = type,
                cluster = cluster_id
            )
        } else {
            econ_multinomial(
                y ~ x1 + x2,
                data = dat,
                vcov = type
            )
        }

        s = summary(fit)
        td = broom::tidy(fit, conf.int = TRUE)

        parameter_names = as.vector(t(outer(
            rownames(s$coefficients),
            colnames(s$coefficients),
            paste,
            sep = ":"
        )))

        tidy_names = paste(td$y.level, td$term, sep = ":")

        idx = match(tidy_names, parameter_names)

        testthat::expect_false(anyNA(idx))
        testthat::expect_equal(length(idx), length(parameter_names))

        testthat::expect_equal(
            td$estimate,
            as.vector(t(s$coefficients))[idx],
            tolerance = 1e-10
        )

        testthat::expect_equal(
            td$std.error,
            as.vector(t(s$standard.errors))[idx],
            tolerance = 1e-10
        )

        testthat::expect_equal(
            td$statistic,
            as.vector(t(s$statistic))[idx],
            tolerance = 1e-10
        )

        testthat::expect_equal(
            td$p.value,
            as.vector(t(s$p.value))[idx],
            tolerance = 1e-10
        )

        critical = if (is.finite(s$df)) {
            stats::qt(0.975, df = s$df)
        } else {
            stats::qnorm(0.975)
        }

        testthat::expect_equal(
            td$conf.low,
            td$estimate - critical * td$std.error,
            tolerance = 1e-10
        )

        testthat::expect_equal(
            td$conf.high,
            td$estimate + critical * td$std.error,
            tolerance = 1e-10
        )
    }
})


testthat::test_that("multinomial tidy validates confidence levels", {

    set.seed(2415)

    dat = data.frame(
        x = rnorm(300),
        y = factor(
            sample(c("A", "B", "C"), 300, replace = TRUE)
        )
    )

    fit = econ_multinomial(y ~ x, data = dat)

    testthat::expect_error(
        broom::tidy(fit, conf.int = TRUE, conf.level = 0),
        "strictly between 0 and 1"
    )

    testthat::expect_error(
        broom::tidy(fit, conf.int = TRUE, conf.level = 1),
        "strictly between 0 and 1"
    )

    testthat::expect_error(
        broom::tidy(fit, conf.int = TRUE, conf.level = NA_real_),
        "strictly between 0 and 1"
    )

    td = broom::tidy(fit)

    testthat::expect_false("conf.low" %in% names(td))
    testthat::expect_false("conf.high" %in% names(td))
})



###### 41. Multinomial Prediction Standard Error Tests ####

testthat::test_that("multinomial augmentation preserves default predictions", {

    set.seed(2417)

    n = 400

    dat = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        region = factor(
            sample(c("East", "North", "South"), n, replace = TRUE)
        ),
        y = factor(
            sample(c("A", "B", "C"), n, replace = TRUE),
            levels = c("A", "B", "C")
        )
    )

    fit = econ_multinomial(
        y ~ x1 + x2 + region,
        data = dat
    )

    aug_default = broom::augment(fit)
    aug_se = broom::augment(fit, se_fit = TRUE)

    prob_columns = paste0(".fitted_", fit$lev)
    se_columns = paste0(".se.fit_", fit$lev)

    testthat::expect_false(
        any(se_columns %in% names(aug_default))
    )

    testthat::expect_true(
        all(se_columns %in% names(aug_se))
    )

    testthat::expect_equal(
        aug_default[, prob_columns],
        aug_se[, prob_columns]
    )

    testthat::expect_identical(
        aug_default$.fitted,
        aug_se$.fitted
    )

    testthat::expect_true(
        all(is.finite(as.matrix(aug_se[, se_columns])))
    )

    testthat::expect_true(
        all(as.matrix(aug_se[, se_columns]) >= 0)
    )

    new_dat = dat[1:10, c("x1", "x2", "region")]

    aug_new = broom::augment(
        fit,
        newdata = new_dat,
        se_fit = TRUE
    )

    testthat::expect_equal(nrow(aug_new), 10L)

    testthat::expect_true(
        all(se_columns %in% names(aug_new))
    )
})


testthat::test_that("multinomial prediction SEs match numerical derivatives", {

    set.seed(2418)

    n = 500

    dat = data.frame(
        x1 = rnorm(n),
        x2 = rnorm(n),
        region = factor(
            sample(c("East", "North", "South"), n, replace = TRUE)
        ),
        cluster_id = sample(seq_len(50), n, replace = TRUE),
        y = factor(
            sample(c("A", "B", "C"), n, replace = TRUE),
            levels = c("A", "B", "C")
        )
    )

    new_dat = dat[1:5, c("x1", "x2", "region")]

    for (type in c("classical", "HC1", "stata")) {

        fit = if (type == "stata") {
            econ_multinomial(
                y ~ x1 + x2 + region,
                data = dat,
                vcov = type,
                cluster = cluster_id
            )
        } else {
            econ_multinomial(
                y ~ x1 + x2 + region,
                data = dat,
                vcov = type
            )
        }

        aug = broom::augment(
            fit,
            newdata = new_dat,
            se_fit = TRUE
        )

        C = coef(fit)

        X = model.matrix(
            delete.response(terms(fit)),
            data = new_dat,
            contrasts.arg = fit$contrasts
        )

        X = X[, colnames(C), drop = FALSE]

        beta = as.vector(t(C))

        parameter_names = unlist(
            lapply(rownames(C), function(r) {
                paste0(r, ":", colnames(C))
            }),
            use.names = FALSE
        )

        V = vcov(fit)[
            parameter_names,
            parameter_names,
            drop = FALSE
        ]

        probability_function = function(b, x_row) {

            coefficients = matrix(
                b,
                nrow = nrow(C),
                byrow = TRUE
            )

            eta = c(
                0,
                as.vector(coefficients %*% x_row)
            )

            eta = eta - max(eta)

            exp(eta) / sum(exp(eta))
        }

        for (i in seq_len(nrow(new_dat))) {

            x_row = as.numeric(X[i, ])

            gradient = matrix(
                0,
                nrow = length(fit$lev),
                ncol = length(beta)
            )

            for (j in seq_along(beta)) {

                h = 1e-5 * max(1, abs(beta[j]))

                beta_plus = beta
                beta_minus = beta

                beta_plus[j] = beta_plus[j] + h
                beta_minus[j] = beta_minus[j] - h

                gradient[, j] = (
                    probability_function(beta_plus, x_row) -
                        probability_function(beta_minus, x_row)
                ) / (2 * h)
            }

            reference_se = sqrt(
                pmax(
                    0,
                    diag(gradient %*% V %*% t(gradient))
                )
            )

            econ_se = as.numeric(
                aug[
                    i,
                    paste0(".se.fit_", fit$lev)
                ]
            )

            testthat::expect_equal(
                econ_se,
                reference_se,
                tolerance = 1e-7
            )
        }
    }
})


testthat::test_that("multinomial augment validates se_fit", {

    set.seed(2419)

    dat = data.frame(
        x = rnorm(300),
        y = factor(
            sample(c("A", "B", "C"), 300, replace = TRUE),
            levels = c("A", "B", "C")
        )
    )

    fit = econ_multinomial(y ~ x, data = dat)

    testthat::expect_error(
        broom::augment(fit, se_fit = NA),
        "'se_fit' must be TRUE or FALSE"
    )

    testthat::expect_error(
        broom::augment(fit, se_fit = 1),
        "'se_fit' must be TRUE or FALSE"
    )

    testthat::expect_error(
        broom::augment(fit, se_fit = c(TRUE, FALSE)),
        "'se_fit' must be TRUE or FALSE"
    )
})



###### 42. Restricted Multinomial Probit Estimation Tests ####

testthat::test_that("multinomial probit probabilities are valid", {

    Omega = matrix(
        c(1, 0.5, 0.5, 1),
        nrow = 2
    )

    eta = matrix(
        0,
        nrow = 5,
        ncol = 3,
        dimnames = list(NULL, c("A", "B", "C"))
    )

    probs = .econ_multinomial_probit_probs(
        eta = eta,
        Omega = Omega
    )

    testthat::expect_equal(
        unname(probs),
        matrix(1 / 3, 5, 3),
        tolerance = 1e-8
    )

    eta = rbind(
        c(0, 1, -1),
        c(0, -1, 1),
        c(0, 0.5, 0.5)
    )

    colnames(eta) = c("A", "B", "C")

    probs = .econ_multinomial_probit_probs(
        eta = eta,
        Omega = Omega
    )

    testthat::expect_equal(
        unname(rowSums(probs)),
        rep(1, nrow(eta)),
        tolerance = 1e-8
    )

    testthat::expect_true(all(probs >= 0))
    testthat::expect_true(all(probs <= 1))
    testthat::expect_gt(probs[1, "B"], probs[1, "C"])
    testthat::expect_gt(probs[2, "C"], probs[2, "B"])
})


testthat::test_that("restricted multinomial probit MLE converges", {

    set.seed(2421)

    n = 120

    dat = data.frame(x = rnorm(n))

    X = model.matrix(~ x, data = dat)

    eta = cbind(
        A = 0,
        B = 0.5 + 0.8 * dat$x,
        C = -0.3 - 0.5 * dat$x
    )

    errors = matrix(
        rnorm(n * 3, sd = sqrt(0.5)),
        nrow = n,
        ncol = 3
    )

    chosen = max.col(
        eta + errors,
        ties.method = "first"
    )

    dat$y = factor(
        c("A", "B", "C")[chosen],
        levels = c("A", "B", "C")
    )

    fit = .econ_multinomial_probit_fit(
        X = X,
        y = dat$y
    )

    testthat::expect_equal(fit$convergence, 0L)

    testthat::expect_equal(
        dim(fit$coefficients),
        c(2L, ncol(X))
    )

    testthat::expect_true(
        all(is.finite(fit$coefficients))
    )

    testthat::expect_equal(
        unname(rowSums(fit$fitted.values)),
        rep(1, n),
        tolerance = 1e-8
    )

    testthat::expect_equal(
        unname(fit$hessian %*% fit$vcov),
        diag(ncol(fit$vcov)),
        tolerance = 1e-8
    )

    testthat::expect_true(
        all(eigen(fit$vcov, symmetric = TRUE)$values > 0)
    )

    observed = fit$fitted.values[
        cbind(seq_len(n), as.integer(dat$y))
    ]

    testthat::expect_equal(
        fit$logLik,
        sum(log(observed)),
        tolerance = 1e-8
    )
})


testthat::test_that("restricted multinomial probit respects weights", {

    set.seed(2422)

    n = 120

    dat = data.frame(
        x = rnorm(n),
        y = factor(
            sample(c("A", "B", "C"), n, replace = TRUE),
            levels = c("A", "B", "C")
        ),
        w = sample(1:3, n, replace = TRUE)
    )

    X = model.matrix(~ x, data = dat)

    fit = .econ_multinomial_probit_fit(
        X = X,
        y = dat$y,
        weights = dat$w
    )

    observed = fit$fitted.values[
        cbind(seq_len(n), as.integer(dat$y))
    ]

    testthat::expect_equal(fit$convergence, 0L)

    testthat::expect_equal(
        unname(fit$weights),
        unname(dat$w)
    )

    testthat::expect_equal(
        fit$logLik,
        sum(dat$w * log(observed)),
        tolerance = 1e-8
    )

    testthat::expect_true(
        all(eigen(fit$vcov, symmetric = TRUE)$values > 0)
    )
})



###### 43. Generalized Multinomial Probit Tests ####

testthat::test_that("generalized multinomial probit probabilities are valid", {

    Omega = matrix(0.5, 3, 3)
    diag(Omega) = 1

    eta = matrix(
        0,
        nrow = 3,
        ncol = 4,
        dimnames = list(NULL, c("A", "B", "C", "D"))
    )

    P = .econ_multinomial_probit_probs_general(
        eta = eta,
        Omega = Omega
    )

    testthat::expect_equal(
        dim(P),
        c(3L, 4L)
    )

    testthat::expect_equal(
        unname(P),
        matrix(0.25, 3, 4),
        tolerance = 1e-8
    )

    testthat::expect_equal(
        unname(rowSums(P)),
        rep(1, 3),
        tolerance = 1e-8
    )

    testthat::expect_true(
        all(P >= 0 & P <= 1)
    )

    eta[1, "B"] = 2

    P_changed = .econ_multinomial_probit_probs_general(
        eta = eta,
        Omega = Omega
    )

    testthat::expect_gt(
        P_changed[1, "B"],
        P[1, "B"]
    )

    testthat::expect_equal(
        unname(rowSums(P_changed)),
        rep(1, 3),
        tolerance = 1e-8
    )
})


testthat::test_that("generalized multinomial probit MLE and covariance are valid", {

    set.seed(1313)

    n = 80
    x = rnorm(n)

    X = model.matrix(~ x)

    eta = cbind(
        A = 0,
        B = 0.4 + 0.5 * x,
        C = -0.3 + 0.2 * x,
        D = 0.1 - 0.4 * x
    )

    errors = matrix(
        rnorm(n * 4, sd = sqrt(0.5)),
        nrow = n,
        ncol = 4
    )

    y = factor(
        c("A", "B", "C", "D")[
            max.col(eta + errors)
        ],
        levels = c("A", "B", "C", "D")
    )

    fit = .econ_multinomial_probit_fit_general(
        X = X,
        y = y
    )

    testthat::expect_equal(
        fit$convergence,
        0L
    )

    testthat::expect_equal(
        dim(fit$coefficients),
        c(3L, 2L)
    )

    testthat::expect_equal(
        dim(fit$vcov),
        c(6L, 6L)
    )

    testthat::expect_true(
        all(is.finite(fit$coefficients))
    )

    testthat::expect_equal(
        unname(rowSums(fit$fitted.values)),
        rep(1, n),
        tolerance = 1e-8
    )

    testthat::expect_true(
        all(eigen(fit$vcov, symmetric = TRUE)$values > 0)
    )

    testthat::expect_equal(
        unname(fit$hessian %*% fit$vcov),
        diag(6),
        tolerance = 1e-8
    )

    observed = fit$fitted.values[
        cbind(seq_len(n), as.integer(y))
    ]

    testthat::expect_equal(
        fit$logLik,
        sum(log(observed)),
        tolerance = 1e-8
    )

    ###### Weighted estimation ####

    set.seed(1314)

    w = sample(1:3, n, replace = TRUE)

    fit_weighted = .econ_multinomial_probit_fit_general(
        X = X,
        y = y,
        weights = w
    )

    testthat::expect_equal(
        fit_weighted$convergence,
        0L
    )

    testthat::expect_equal(
        unname(fit_weighted$weights),
        unname(w)
    )

    observed_weighted = fit_weighted$fitted.values[
        cbind(seq_len(n), as.integer(y))
    ]

    testthat::expect_equal(
        fit_weighted$logLik,
        sum(w * log(observed_weighted)),
        tolerance = 1e-8
    )
})



###### 44. Multinomial Probit Robust Covariance Tests ####

testthat::test_that("probit scores and HC covariance are valid", {

    ###### 1. Generate four-category data ####

    set.seed(1313)

    n = 80
    x = rnorm(n)

    X = model.matrix(~ x)

    eta = cbind(
        A = 0,
        B = 0.4 + 0.5 * x,
        C = -0.3 + 0.2 * x,
        D = 0.1 - 0.4 * x
    )

    errors = matrix(
        rnorm(n * 4, sd = sqrt(0.5)),
        nrow = n,
        ncol = 4
    )

    y = factor(
        c("A", "B", "C", "D")[
            max.col(eta + errors)
        ],
        levels = c("A", "B", "C", "D")
    )

    fit = .econ_multinomial_probit_fit_general(
        X = X,
        y = y
    )

    ###### 2. Validate individual scores ####

    S = .econ_multinomial_probit_scores(
        fit = fit,
        X = X,
        y = y
    )

    testthat::expect_equal(
        dim(S),
        c(n, 6L)
    )

    testthat::expect_identical(
        colnames(S),
        colnames(fit$vcov)
    )

    testthat::expect_true(
        all(is.finite(S))
    )

    testthat::expect_lt(
        max(abs(colSums(S))),
        0.1
    )

    ###### 3. Validate classical and HC covariance ####

    V = fit$vcov
    k = ncol(S)

    V_classical = .econ_multinomial_probit_vcov(
        fit, X, y,
        vcov = "classical"
    )

    V_HC0 = .econ_multinomial_probit_vcov(
    fit = fit,
    X = X,
    y = y,
    vcov = "HC0",
    scores = S
)

V_HC1 = .econ_multinomial_probit_vcov(
    fit = fit,
    X = X,
    y = y,
    vcov = "HC1",
    scores = S
)

    HC0_reference = V %*% crossprod(S) %*% V
    HC1_reference = HC0_reference * n / (n - k)

    testthat::expect_equal(
        V_classical,
        V
    )

    testthat::expect_equal(
        unname(V_HC0),
        unname(HC0_reference),
        tolerance = 1e-8
    )

    testthat::expect_equal(
        unname(V_HC1),
        unname(HC1_reference),
        tolerance = 1e-8
    )

    ###### 4. Validate clustered covariance ####

    cluster = rep(seq_len(20), length.out = n)

    G = length(unique(cluster))

    U = rowsum(
        S,
        group = cluster,
        reorder = TRUE
    )

    CR0_reference = V %*% crossprod(U) %*% V

    corrections = c(
        CR0 = 1,
        CR1 = (n - 1) / (n - k),
        CR1p = G / (G - k),
        CR1S = G * (n - 1) / ((G - 1) * (n - k)),
        stata = G * (n - 1) / ((G - 1) * (n - k))
    )

    for (type in names(corrections)) {

    result = .econ_multinomial_probit_vcov(
        fit = fit,
        X = X,
        y = y,
        vcov = type,
        cluster = cluster,
        scores = S
    )

    reference = CR0_reference * corrections[[type]]

    testthat::expect_equal(
        unname(result),
        unname(reference),
        tolerance = 1e-8,
        info = paste("Covariance type:", type)
    )
}

    ###### 5. Validate weighted covariance ####

    set.seed(1314)

    w = sample(1:3, n, replace = TRUE)

    fit_weighted = .econ_multinomial_probit_fit_general(
        X = X,
        y = y,
        weights = w
    )

    S_weighted = .econ_multinomial_probit_scores(
        fit = fit_weighted,
        X = X,
        y = y
    )

    V_weighted = .econ_multinomial_probit_vcov(
    fit = fit_weighted,
    X = X,
    y = y,
    vcov = "HC0",
    scores = S_weighted
)

    weighted_reference =
        fit_weighted$vcov %*%
        crossprod(S_weighted) %*%
        fit_weighted$vcov

    testthat::expect_true(
        all(is.finite(S_weighted))
    )

    testthat::expect_equal(
        unname(V_weighted),
        unname(weighted_reference),
        tolerance = 1e-8
    )

    ###### 6. Validate cluster error handling ####

    testthat::expect_error(
        .econ_multinomial_probit_vcov(
            fit = fit,
            X = X,
            y = y,
            vcov = "CR0"
        ),
        "cluster"
    )

    testthat::expect_error(
        .econ_multinomial_probit_vcov(
            fit = fit,
            X = X,
            y = y,
            vcov = "CR1p",
            cluster = rep(1:3, length.out = n)
        ),
        "G > k"
    )

###### 7. Validate Direct and Reused Scores ####

V_direct = .econ_multinomial_probit_vcov(
    fit = fit,
    X = X,
    y = y,
    vcov = "HC0"
)

testthat::expect_equal(
    V_direct,
    V_HC0,
    tolerance = 1e-10
)

})



###### 45. Public Multinomial Probit API Tests ####

testthat::test_that("public multinomial probit estimation and prediction work", {

    ###### 1. Generate test data ####

    set.seed(1515)

    n = 72L
    x = rnorm(n)

    eta = cbind(
        A = 0,
        B = 0.5 + 0.7 * x,
        C = -0.2 - 0.4 * x
    )

    errors = matrix(
        rnorm(n * 3L, sd = sqrt(0.5)),
        nrow = n,
        ncol = 3L
    )

    dat = data.frame(
        x = x,
        y = factor(
            c("A", "B", "C")[
                max.col(eta + errors)
            ],
            levels = c("A", "B", "C")
        ),
        w = sample(1:3, n, replace = TRUE),
        cluster = rep(seq_len(12L), each = 6L)
    )

    ###### 2. Public estimation ####

    fit = econ_multinomial(
        y ~ x,
        data = dat,
        link = "probit"
    )

    testthat::expect_s3_class(
        fit,
        "econ_multinomial"
    )

    testthat::expect_identical(
        fit$link,
        "probit"
    )

    testthat::expect_false(
        inherits(fit, "multinom")
    )

    testthat::expect_equal(
        dim(stats::coef(fit)),
        c(2L, 2L)
    )

    testthat::expect_equal(
        dim(stats::vcov(fit)),
        c(4L, 4L)
    )

    testthat::expect_true(
        is.finite(as.numeric(stats::logLik(fit)))
    )

    testthat::expect_equal(
        stats::nobs(fit),
        n
    )

    ###### 3. Predictions ####

    P = stats::predict(fit, type = "probs")

    testthat::expect_equal(
        dim(P),
        c(n, 3L)
    )

    testthat::expect_equal(
        unname(rowSums(P)),
        rep(1, n),
        tolerance = 1e-8
    )

    testthat::expect_true(
        all(is.finite(P) & P >= 0 & P <= 1)
    )

    testthat::expect_equal(
        stats::fitted(fit),
        P
    )

    predicted_class = stats::predict(
        fit,
        type = "class"
    )

    testthat::expect_identical(
        as.character(predicted_class),
        colnames(P)[max.col(P, ties.method = "first")]
    )

    newdata = data.frame(
        x = c(-1, 0, 1)
    )

    P_new = stats::predict(
        fit,
        newdata = newdata,
        type = "probs"
    )

    testthat::expect_equal(
        dim(P_new),
        c(3L, 3L)
    )

    testthat::expect_equal(
        unname(rowSums(P_new)),
        rep(1, 3),
        tolerance = 1e-8
    )

    testthat::expect_identical(
        colnames(P_new),
        levels(dat$y)
    )

    ###### 4. Scores ####

    S = estfun.econ_multinomial(fit)

    testthat::expect_equal(
        dim(S),
        c(n, 4L)
    )

    testthat::expect_true(
        all(is.finite(S))
    )

    testthat::expect_identical(
        colnames(S),
        colnames(stats::vcov(fit))
    )

    ###### 5. Classical covariance ####

    testthat::expect_equal(
        stats::vcov(fit),
        fit$vcov,
        tolerance = 1e-10
    )

    ###### 6. HC1 covariance ####

    fit_hc = fit
    fit_hc$vcov_type = "HC1"

    V_hc = stats::vcov(fit_hc)

    testthat::expect_equal(
        dim(V_hc),
        c(4L, 4L)
    )

    testthat::expect_true(
        all(is.finite(V_hc))
    )

    testthat::expect_identical(
        fit_hc$vcov_type,
        "HC1"
    )

    ###### 7. Cluster-robust covariance ####

    fit_cr = fit
    fit_cr$vcov_type = "CR1S"
    fit_cr$cluster = dat$cluster

    V_cr = stats::vcov(fit_cr)

    testthat::expect_equal(
        dim(V_cr),
        c(4L, 4L)
    )

    testthat::expect_true(
        all(is.finite(V_cr))
    )

    testthat::expect_identical(
        fit_cr$cluster,
        dat$cluster
    )

    ###### 8. Summary and clustered inference ####

    s = summary(fit_cr)

    testthat::expect_s3_class(
        s,
        "summary.econ_multinomial"
    )

    testthat::expect_identical(
        s$link,
        "probit"
    )

    testthat::expect_equal(
        s$df,
        11
    )

    testthat::expect_equal(
        as.vector(t(s$standard.errors)),
        unname(sqrt(diag(V_cr))),
        ignore_attr = TRUE,
        tolerance = 1e-10
    )

    expected_p = 2 * stats::pt(
        -abs(s$statistic),
        df = 11
    )

    testthat::expect_equal(
        s$p.value,
        expected_p,
        tolerance = 1e-10
    )

    ###### 9. Confidence intervals ####

    CI = stats::confint(fit_cr)

    testthat::expect_equal(
        dim(CI),
        c(4L, 2L)
    )

    testthat::expect_true(
        all(is.finite(CI))
    )

    testthat::expect_identical(
        rownames(CI),
        colnames(V_cr)
    )

    selected = stats::confint(
        fit_cr,
        parm = "B:x"
    )

    testthat::expect_equal(
        dim(selected),
        c(1L, 2L)
    )

    testthat::expect_identical(
        rownames(selected),
        "B:x"
    )

    ###### 10. Broom reporting ####

    tidy_result = broom::tidy(fit_cr)

    testthat::expect_equal(
        nrow(tidy_result),
        4L
    )

    testthat::expect_true(
        all(is.finite(tidy_result$std.error))
    )

    glance_result = broom::glance(fit)

    testthat::expect_true(
        is.finite(glance_result$pseudo.r.squared)
    )

    ###### 11. Preserve multinomial logit ####

    fit_logit = econ_multinomial(
        y ~ x,
        data = dat,
        link = "logit"
    )

    testthat::expect_true(
        inherits(fit_logit, "multinom")
    )

    testthat::expect_identical(
        fit_logit$link,
        "logit"
    )

    testthat::expect_equal(
        unname(rowSums(
            stats::predict(fit_logit, type = "probs")
        )),
        rep(1, n),
        tolerance = 1e-8
    )
})



###### Shared Fixed-Probit Fixture for Sections 46 and 47 ####

.stage46_47_fixture = local({

    cached = NULL

    function() {

        if (!is.null(cached)) return(cached)

        set.seed(1616)
        n = 60L

        dat = data.frame(
            x = rnorm(n),
            group = factor(rep(c("Control", "Treatment"), length.out = n)),
            w = sample(1:3, n, replace = TRUE),
            cl = rep(seq_len(15L), each = 4L)
        )

        eta = cbind(
            A = 0,
            B = 0.3 + 0.7 * dat$x + 0.3 * (dat$group == "Treatment"),
            C = -0.2 - 0.5 * dat$x + 0.4 * (dat$group == "Treatment")
        )

        errors = matrix(
            rnorm(n * 3L, sd = sqrt(0.5)),
            nrow = n,
            ncol = 3L
        )

        dat$y = factor(
            c("A", "B", "C")[max.col(eta + errors)],
            levels = c("A", "B", "C")
        )

        fit = econ_multinomial(
            y ~ x + group,
            data = dat,
            link = "probit"
        )

        fit_w = econ_multinomial(
            y ~ x + group,
            data = dat,
            weights = w,
            link = "probit",
            vcov = "CR1S",
            cluster = cl
        )

        cached <<- list(dat = dat, fit = fit, fit_w = fit_w)
        cached
    }
})


###### 46. Multinomial Probit Prediction SE and Marginal Effects Tests ####

testthat::test_that("probit prediction SEs and augmentation work", {

    ###### 1. Reuse reproducible fixed-probit fit ####

    z = .stage46_47_fixture()
    dat = z$dat
    fit = z$fit
    n = nrow(dat)

    ###### 3. In-sample augmentation ####

    aug = broom::augment(
        fit,
        se_fit = TRUE
    )

    outcomes = levels(dat$y)

    se_columns = paste0(
        ".se.fit_",
        outcomes
    )

    probability_columns = paste0(
        ".fitted_",
        outcomes
    )

    testthat::expect_equal(
        nrow(aug),
        n
    )

    testthat::expect_true(
        all(se_columns %in% names(aug))
    )

    testthat::expect_true(
        all(probability_columns %in% names(aug))
    )

    SE = as.matrix(
        aug[, se_columns, drop = FALSE]
    )

    P = as.matrix(
        aug[, probability_columns, drop = FALSE]
    )

    testthat::expect_true(
        all(is.finite(SE))
    )

    testthat::expect_true(
        all(SE >= 0)
    )

    testthat::expect_equal(
        unname(rowSums(P)),
        rep(1, n),
        tolerance = 1e-8
    )

    testthat::expect_equal(
        unname(P),
        unname(predict(fit, type = "probs")),
        tolerance = 1e-8
    )

    ###### 4. Out-of-sample augmentation ####

    newdata = data.frame(
        x = c(-1, 0, 1),
        group = factor(
            c("Control", "Treatment", "Control"),
            levels = levels(dat$group)
        )
    )

    aug_new = broom::augment(
        fit,
        newdata = newdata,
        se_fit = TRUE
    )

    testthat::expect_equal(
        nrow(aug_new),
        3L
    )

    testthat::expect_true(
        all(se_columns %in% names(aug_new))
    )

    testthat::expect_true(
        all(is.finite(
            as.matrix(aug_new[, se_columns, drop = FALSE])
        ))
    )

    testthat::expect_true(
        all(
            as.matrix(aug_new[, se_columns, drop = FALSE]) >= 0
        )
    )

    testthat::expect_equal(
        unname(as.matrix(
            aug_new[, probability_columns, drop = FALSE]
        )),
        unname(predict(
            fit,
            newdata = newdata,
            type = "probs"
        )),
        tolerance = 1e-8
    )
})


testthat::test_that("probit continuous and categorical AMEs work", {

    ###### 1. Generate reproducible data ####

    ###### 1. Reuse reproducible fixed-probit fit ####

    z = .stage46_47_fixture()
    dat = z$dat
    fit = z$fit
    n = nrow(dat)

    ###### 3. Continuous marginal effects ####

    ame_x = marginal_effects.econ_multinomial(
        fit,
        variables = "x"
    )

    testthat::expect_identical(
        as.character(ame_x$outcome),
        levels(dat$y)
    )

    testthat::expect_true(
        all(is.finite(ame_x$AME))
    )

    testthat::expect_lt(
        abs(sum(ame_x$AME)),
        1e-5
    )

    testthat::expect_true(
        all(is.finite(ame_x$std.error))
    )

    testthat::expect_true(
        all(ame_x$std.error >= 0)
    )

    ###### 4. Categorical discrete probability changes ####

    ame_group = marginal_effects.econ_multinomial(
        fit,
        variables = "group"
    )

    testthat::expect_identical(
        as.character(ame_group$outcome),
        levels(dat$y)
    )

    testthat::expect_true(
        all(is.finite(ame_group$AME))
    )

    testthat::expect_lt(
        abs(sum(ame_group$AME)),
        1e-5
    )

    testthat::expect_true(
        all(
            ame_group$contrast == "Treatment - Control"
        )
    )

    ###### 5. Weighted and clustered marginal effects ####

    fit_weighted = z$fit_w

    ame_weighted = marginal_effects.econ_multinomial(
        fit_weighted,
        variables = c("x", "group"),
        weighted = TRUE
    )

    testthat::expect_equal(
        nrow(ame_weighted),
        6L
    )

    testthat::expect_true(
        all(is.finite(ame_weighted$AME))
    )

    testthat::expect_true(
        all(is.finite(ame_weighted$std.error))
    )

    testthat::expect_true(
        all(ame_weighted$std.error >= 0)
    )

    testthat::expect_true(
        all(is.finite(ame_weighted$p.value))
    )

    testthat::expect_true(
        all(
            ame_weighted$p.value >= 0 &
                ame_weighted$p.value <= 1
        )
    )

    ###### 6. Preserve multinomial logit behavior ####

    fit_logit = econ_multinomial(
        y ~ x + group,
        data = dat,
        link = "logit"
    )

    aug_logit = broom::augment(
        fit_logit,
        se_fit = TRUE
    )

    ame_logit = marginal_effects.econ_multinomial(
        fit_logit,
        variables = "x"
    )

    testthat::expect_equal(
        nrow(aug_logit),
        n
    )

    testthat::expect_true(
        all(is.finite(as.matrix(
            aug_logit[
                ,
                paste0(".se.fit_", levels(dat$y)),
                drop = FALSE
            ]
        )))
    )

    testthat::expect_true(
        all(is.finite(ame_logit$AME))
    )
})



###### 47. Multinomial Probit Numerical Accuracy Tests ####

testthat::test_that("probit prediction SEs and AMEs match numerical references", {

    ###### 1. Generate data ####

    z = .stage46_47_fixture()
    dat = z$dat
    fit = z$fit
    fit_w = z$fit_w
    n = nrow(dat)

    ###### 2. Independent prediction function ####

    predict_beta = function(model, beta, data) {

        modified = model

        modified$coefficients = matrix(
            beta,
            nrow = nrow(coef(model)),
            ncol = ncol(coef(model)),
            byrow = TRUE,
            dimnames = dimnames(coef(model))
        )

        as.matrix(predict(
            modified,
            newdata = data,
            type = "probs"
        ))
    }

    ###### 3. Prediction standard errors ####

    newdata = data.frame(
        x = c(-1, 0, 1),
        group = factor(
            c("Control", "Treatment", "Control"),
            levels = levels(dat$group)
        )
    )

    beta = as.vector(t(coef(fit)))
    V = vcov(fit)

    P = predict_beta(fit, beta, newdata)

    gradients = array(
        0,
        dim = c(nrow(newdata), ncol(P), length(beta))
    )

    for (k in seq_along(beta)) {

        h = 1e-5 * max(1, abs(beta[k]))

        plus = beta
        minus = beta

        plus[k] = plus[k] + h
        minus[k] = minus[k] - h

        gradients[, , k] = (
            predict_beta(fit, plus, newdata) -
            predict_beta(fit, minus, newdata)
        ) / (2 * h)
    }

    SE_reference = matrix(
        0,
        nrow = nrow(newdata),
        ncol = ncol(P)
    )

    for (j in seq_len(ncol(P))) {

        G = matrix(
            gradients[, j, ],
            nrow = nrow(newdata),
            ncol = length(beta)
        )

        SE_reference[, j] = sqrt(pmax(
            0,
            rowSums((G %*% V) * G)
        ))
    }

    aug = broom::augment(
        fit,
        newdata = newdata,
        se_fit = TRUE
    )

    SE_econ = as.matrix(
        aug[, paste0(".se.fit_", colnames(P)), drop = FALSE]
    )

    testthat::expect_equal(
        unname(SE_econ),
        unname(SE_reference),
        tolerance = 1e-5
    )

    ###### 4. Independent continuous and factor effects ####

    continuous_ame = function(model, beta, data, weights) {

        h = 1e-4

        plus = data
        minus = data

        plus$x = plus$x + h
        minus$x = minus$x - h

        effect = (
            predict_beta(model, beta, plus) -
            predict_beta(model, beta, minus)
        ) / (2 * h)

        colSums(effect * (weights / sum(weights)))
    }

    factor_ame = function(model, beta, data, weights) {

        high = data
        low = data

        high$group = factor(
            rep("Treatment", nrow(data)),
            levels = levels(data$group)
        )

        low$group = factor(
            rep("Control", nrow(data)),
            levels = levels(data$group)
        )

        effect = predict_beta(model, beta, high) -
            predict_beta(model, beta, low)

        colSums(effect * (weights / sum(weights)))
    }

    equal_weights = rep(1, n)

    reference_x = continuous_ame(
        fit, beta, dat, equal_weights
    )

    reference_group = factor_ame(
        fit, beta, dat, equal_weights
    )

    ame_x = marginal_effects.econ_multinomial(
        fit,
        variables = "x"
    )

    ame_group = marginal_effects.econ_multinomial(
        fit,
        variables = "group"
    )

    testthat::expect_equal(
        unname(ame_x$AME),
        unname(reference_x),
        tolerance = 1e-5
    )

    testthat::expect_equal(
        unname(ame_group$AME),
        unname(reference_group),
        tolerance = 1e-6
    )

    ###### 5. Independently weighted marginal effects ####

    beta_w = as.vector(t(coef(fit_w)))

    weighted_reference = c(
        continuous_ame(fit_w, beta_w, dat, dat$w),
        factor_ame(fit_w, beta_w, dat, dat$w)
    )

    ame_w = marginal_effects.econ_multinomial(
        fit_w,
        variables = c("x", "group"),
        weighted = TRUE
    )

    testthat::expect_equal(
        unname(ame_w$AME),
        unname(weighted_reference),
        tolerance = 1e-5
    )

    ###### 6. Independently calculated delta-method SEs ####

    ame_function = function(b) {

        c(
            continuous_ame(fit_w, b, dat, dat$w),
            factor_ame(fit_w, b, dat, dat$w)
        )
    }

    Jacobian = matrix(
        0,
        nrow = length(weighted_reference),
        ncol = length(beta_w)
    )

    for (k in seq_along(beta_w)) {

        h = 1e-4 * max(1, abs(beta_w[k]))

        plus = beta_w
        minus = beta_w

        plus[k] = plus[k] + h
        minus[k] = minus[k] - h

        Jacobian[, k] = (
            ame_function(plus) -
            ame_function(minus)
        ) / (2 * h)
    }

    V_w = vcov(fit_w)

    SE_reference_ame = sqrt(pmax(
        0,
        diag(Jacobian %*% V_w %*% t(Jacobian))
    ))

    testthat::expect_equal(
        unname(ame_w$std.error),
        unname(SE_reference_ame),
        tolerance = 1e-4
    )

    ###### 7. Probability identities ####

    testthat::expect_lt(
        abs(sum(reference_x)),
        1e-6
    )

    testthat::expect_lt(
        abs(sum(reference_group)),
        1e-6
    )

    testthat::expect_lt(
        abs(sum(tail(weighted_reference, 3L))),
        1e-6
    )
})


###### 48. Estimable Multinomial Probit Error Covariance ####


###### 48.1. Identified Cholesky Covariance ####

testthat::test_that("Estimated probit covariance is identified and positive definite", {

    for (M in 2:5) {

        q = M * (M + 1L) / 2L - 1L

        gamma = rep(0, q)

        Omega = .econ_multinomial_probit_omega(
            gamma = gamma,
            M = M
        )

        testthat::expect_equal(
            dim(Omega),
            c(M, M)
        )

        testthat::expect_equal(
            Omega,
            t(Omega),
            tolerance = 1e-12
        )

        testthat::expect_equal(
            Omega[1L, 1L],
            1,
            tolerance = 1e-12
        )

        testthat::expect_true(
            all(eigen(Omega, symmetric = TRUE)$values > 0)
        )

        testthat::expect_length(
            .econ_multinomial_probit_cov_names(M),
            q
        )

        testthat::expect_equal(
            length(unique(.econ_multinomial_probit_cov_names(M))),
            q
        )
    }

    gamma = c(
        0.5,
        log(sqrt(0.75))
    )

    testthat::expect_equal(
        unname(.econ_multinomial_probit_omega(gamma, 2L)),
        matrix(c(1, 0.5, 0.5, 1), 2L),
        tolerance = 1e-12
    )

    testthat::expect_error(
        .econ_multinomial_probit_omega(0, 2L),
        "Invalid"
    )

    testthat::expect_error(
        .econ_multinomial_probit_omega(c(0, Inf), 2L),
        "Invalid"
    )
})


###### 48.2. Independent Probability Identities ####

testthat::test_that("Estimated probit probabilities agree with direct MVN integrals", {

    eta = rbind(
        c(0, 0.4, -0.3),
        c(0, -0.2, 0.7),
        c(0, 0, 0)
    )

    colnames(eta) = c("A", "B", "C")

    Omega = matrix(
        c(1, 0.2, 0.2, 1.4),
        2L
    )

    P = .econ_multinomial_probit_probs_general(
        eta = eta,
        Omega = Omega
    )

    Sigma = matrix(0, 3L, 3L)
    Sigma[2:3, 2:3] = Omega

    reference = matrix(
        NA_real_,
        3L,
        3L
    )

    for (j in 1:3) {

        others = setdiff(1:3, j)

        D = diag(3L)[others, , drop = FALSE] -
            matrix(
                diag(3L)[j, ],
                nrow = 2L,
                ncol = 3L,
                byrow = TRUE
            )

        S = D %*% Sigma %*% t(D)

        for (i in 1:3) {

            reference[i, j] = as.numeric(
                mvtnorm::pmvnorm(
                    upper = eta[i, j] - eta[i, others],
                    sigma = S,
                    algorithm = mvtnorm::Miwa()
                )
            )
        }
    }

    testthat::expect_equal(
        unname(P),
        reference,
        tolerance = 1e-9
    )

    testthat::expect_equal(
        rowSums(P),
        rep(1, nrow(P)),
        tolerance = 1e-8
    )

    testthat::expect_true(
        all(P > 0 & P < 1)
    )
})


###### 48.3. Shared Simulation Data ####

.stage18_fixture = local({

    cached = NULL

    function() {

        if (!is.null(cached)) {
            return(cached)
        }

        set.seed(1818)

        n = 300L

        dat = data.frame(
            x = rnorm(n),
            cl = rep(
                seq_len(n / 10L),
                each = 10L
            )
        )

        z1 = rnorm(n)
        z2 = rnorm(n)

        error_B = z1

        error_C = 0.35 * z1 +
            sqrt(1.25) * z2

        utilities = cbind(
            A = 0,
            B = 0.25 + 0.8 * dat$x + error_B,
            C = -0.2 - 0.5 * dat$x + error_C
        )

        dat$y = factor(
            colnames(utilities)[max.col(utilities)],
            levels = c("A", "B", "C")
        )

        fit = econ_multinomial(
            y ~ x,
            data = dat,
            link = "probit",
            error_covariance = "estimated"
        )

        cached <<- list(
            dat = dat,
            fit = fit
        )

        cached
    }
})


###### 48.4. Joint Parameter Estimation ####

testthat::test_that("Estimated probit stores identified covariance parameters", {

    z = .stage18_fixture()
    fit = z$fit

    testthat::expect_identical(
        fit$error_covariance,
        "estimated"
    )

    testthat::expect_equal(
        fit$convergence,
        0L
    )

    testthat::expect_length(
        fit$covariance_parameters,
        2L
    )

    testthat::expect_equal(
        dim(fit$Omega),
        c(2L, 2L)
    )

    testthat::expect_equal(
        fit$Omega[1L, 1L],
        1,
        tolerance = 1e-10
    )

    testthat::expect_true(
        all(eigen(fit$Omega, symmetric = TRUE)$values > 0)
    )

    testthat::expect_equal(
        unname(fit$Omega),
        unname(
            .econ_multinomial_probit_omega(
                fit$covariance_parameters,
                2L
            )
        ),
        tolerance = 1e-10
    )

    testthat::expect_equal(
        ncol(vcov(fit)),
        length(coef(fit)) + 2L
    )

    testthat::expect_identical(
        colnames(vcov(fit)),
        fit$parameter_names
    )

    testthat::expect_equal(
        unname(fit$hessian %*% fit$vcov),
        diag(ncol(fit$vcov)),
        tolerance = 1e-5
    )
})


###### 48.5. Log-Likelihood and Information Criteria ####

testthat::test_that("Estimated covariance parameters count toward information criteria", {

    z = .stage18_fixture()
    fit = z$fit

    y = stats::model.response(
        fit$econ_model
    )

    p = predict(
        fit,
        type = "probs"
    )

    idx = cbind(
        seq_len(nrow(p)),
        as.integer(y)
    )

    ll_reference = sum(
        fit$econ_weights *
            log(pmax(p[idx], .Machine$double.xmin))
    )

    k = length(coef(fit)) +
        length(fit$covariance_parameters)

    testthat::expect_equal(
        as.numeric(logLik(fit)),
        ll_reference,
        tolerance = 1e-6
    )

    testthat::expect_equal(
        attr(logLik(fit), "df"),
        k
    )

    testthat::expect_equal(
        AIC(fit),
        -2 * ll_reference + 2 * k,
        tolerance = 1e-6
    )

    testthat::expect_equal(
        BIC(fit),
        -2 * ll_reference + log(nobs(fit)) * k,
        tolerance = 1e-6
    )

    testthat::expect_equal(
        as.numeric(glance.econ_multinomial(fit)$edf),
        k
    )
})


###### 48.6. Scores and Robust Covariance ####

testthat::test_that("Joint probit scores produce HC and clustered covariance", {

    z = .stage18_fixture()
    fit = z$fit

    X = stats::model.matrix(
        fit$terms,
        fit$econ_model,
        contrasts.arg = fit$contrasts
    )

    y = stats::model.response(
        fit$econ_model
    )

    S = .econ_multinomial_probit_scores(
        fit,
        X,
        y
    )

    V = fit$vcov

    testthat::expect_equal(
        dim(S),
        c(nrow(X), ncol(V))
    )

    testthat::expect_identical(
        colnames(S),
        colnames(V)
    )

    testthat::expect_true(
        all(is.finite(S))
    )

    hc0_reference = V %*%
        crossprod(S) %*%
        V

    hc0 = .econ_multinomial_probit_vcov(
        fit,
        X,
        y,
        vcov = "HC0",
        scores = S
    )

    hc1 = .econ_multinomial_probit_vcov(
        fit,
        X,
        y,
        vcov = "HC1",
        scores = S
    )

    testthat::expect_equal(
        unname(hc0),
        unname(hc0_reference),
        tolerance = 1e-7
    )

    testthat::expect_equal(
        unname(hc1),
        unname(
            hc0 * nrow(S) /
                (nrow(S) - ncol(S))
        ),
        tolerance = 1e-7
    )

    cl = z$dat$cl

    U = rowsum(
        S,
        cl
    )

    cr0_reference = V %*%
        crossprod(U) %*%
        V

    cr0 = .econ_multinomial_probit_vcov(
        fit,
        X,
        y,
        vcov = "CR0",
        cluster = cl,
        scores = S
    )

    testthat::expect_equal(
        unname(cr0),
        unname(cr0_reference),
        tolerance = 1e-7
    )

    G = length(unique(cl))

    cr1s = .econ_multinomial_probit_vcov(
        fit,
        X,
        y,
        vcov = "CR1S",
        cluster = cl,
        scores = S
    )

    testthat::expect_equal(
        unname(cr1s),
        unname(
            cr0 * G * (nrow(S) - 1) /
                ((G - 1) * (nrow(S) - ncol(S)))
        ),
        tolerance = 1e-7
    )
})


###### 48.7. Inference and Confidence Intervals ####

testthat::test_that("Probit coefficient inference uses joint parameter uncertainty", {

    fit = .stage18_fixture()$fit

    s = summary(fit)

    V = vcov(fit)

    C = coef(fit)

    beta_names = unlist(
        lapply(
            rownames(C),
            function(v) {
                paste0(v, ":", colnames(C))
            }
        ),
        use.names = FALSE
    )

    testthat::expect_equal(
        unname(
            as.vector(t(s$standard.errors))
        ),
        unname(
            sqrt(diag(
                V[
                    beta_names,
                    beta_names,
                    drop = FALSE
                ]
            ))
        ),
        tolerance = 1e-8
    )

    ci = confint(fit)

    testthat::expect_identical(
        rownames(ci),
        beta_names
    )

    testthat::expect_equal(
        nrow(tidy.econ_multinomial(fit)),
        length(C)
    )

    testthat::expect_equal(
        dim(C),
        dim(s$coefficients)
    )
})


###### 48.8. Predictions and Standard Errors ####

testthat::test_that("Predictions use estimated covariance and joint delta-method inference", {

    fit = .stage18_fixture()$fit

    newdata = data.frame(
        x = c(-0.8, 0, 0.8)
    )

    X = model.matrix(
        stats::delete.response(fit$terms),
        newdata,
        contrasts.arg = fit$contrasts
    )

    eta = cbind(
        0,
        X %*% t(coef(fit))
    )

    colnames(eta) = fit$lev

    p_reference = .econ_multinomial_probit_probs_general(
        eta = eta,
        Omega = fit$Omega,
        normalize = TRUE
    )

    p = predict(
        fit,
        newdata = newdata,
        type = "probs"
    )

    testthat::expect_equal(
        unname(p),
        unname(p_reference),
        tolerance = 1e-9
    )

    testthat::expect_equal(
        rowSums(p),
        rep(1, nrow(p)),
        tolerance = 1e-8
    )

    se = .econ_multinomial_probit_prediction_se(
        fit,
        newdata
    )

    testthat::expect_equal(
        dim(se),
        dim(p)
    )

    testthat::expect_true(
        all(is.finite(se) & se >= 0)
    )

    aug = augment.econ_multinomial(
        fit,
        newdata = newdata,
        se_fit = TRUE
    )

    for (j in seq_along(fit$lev)) {

        nm = paste0(
            ".se.fit_",
            make.names(fit$lev[j])
        )

        testthat::expect_equal(
            unname(aug[[nm]]),
            unname(se[, j]),
            tolerance = 1e-6
        )
    }
})


###### 48.9. Average Marginal Effects ####

testthat::test_that("Estimated covariance probit produces finite AMEs", {

    fit = .stage18_fixture()$fit

    ame = marginal_effects.econ_multinomial(
        fit,
        variables = "x"
    )

    testthat::expect_equal(
        nrow(ame),
        length(fit$lev)
    )

    testthat::expect_true(
        all(is.finite(ame$AME))
    )

    testthat::expect_true(
        all(
            is.finite(ame$std.error) &
                ame$std.error >= 0
        )
    )

    testthat::expect_lt(
        abs(sum(ame$AME)),
        1e-5
    )
})


###### 48.10. Multinomial Logit Compatibility ####

testthat::test_that("Estimated covariance is rejected for multinomial logit", {

    dat = data.frame(
        y = factor(
            rep(c("A", "B", "C"), each = 10)
        ),
        x = rnorm(30)
    )

    testthat::expect_error(
        econ_multinomial(
            y ~ x,
            data = dat,
            link = "logit",
            error_covariance = "estimated"
        ),
        "only supported"
    )
})


###### 49. Estimated-Covariance Probit Sample Handling ####


###### 49.1. Shared Weighted and Incomplete Sample ####

.stage49_fixture = local({

    cached = NULL

    function() {

        if (!is.null(cached)) {
            return(cached)
        }

        dat = .stage18_fixture()$dat

        set.seed(1949)

        dat$w = runif(
            nrow(dat),
            min = 0.75,
            max = 1.25
        )

        dat$keep = rep(
            c(TRUE, TRUE, TRUE, TRUE, FALSE),
            length.out = nrow(dat)
        )

        dat$x[c(5, 20, 45, 100, 155)] = NA_real_

        expected = which(
            dat$keep &
                stats::complete.cases(
                    dat[, c("y", "x", "w")]
                )
        )

        fit = econ_multinomial(
            y ~ x,
            data = dat,
            weights = w,
            subset = keep,
            link = "probit",
            error_covariance = "estimated",
            vcov = "classical",
            cluster = cl
        )

        cached <<- list(
            dat = dat,
            fit = fit,
            expected = expected
        )

        cached
    }
})


###### 49.2. Missing Covariates and Subset Alignment ####

testthat::test_that(
    "Estimated probit aligns observations after missingness and subsetting",
    {

        z = .stage49_fixture()

        fit = z$fit
        dat = z$dat
        expected = z$expected

        testthat::expect_equal(
            fit$used,
            expected
        )

        testthat::expect_equal(
            fit$omitted,
            setdiff(seq_len(nrow(dat)), expected)
        )

        testthat::expect_equal(
            nobs(fit),
            length(expected)
        )

        testthat::expect_equal(
            nrow(fit$econ_model),
            length(expected)
        )

        testthat::expect_equal(
            nrow(fitted(fit)),
            length(expected)
        )

        testthat::expect_true(
            all(is.finite(fitted(fit)))
        )

        testthat::expect_equal(
            rowSums(fitted(fit)),
            rep(1, length(expected)),
            tolerance = 1e-7
        )
    }
)


###### 49.3. Observation Weights and Weighted Likelihood ####

testthat::test_that(
    "Estimated probit preserves weights and the weighted likelihood",
    {

        z = .stage49_fixture()

        fit = z$fit
        dat = z$dat
        expected = z$expected

        ###### 1. Validate observation weights ####

        testthat::expect_equal(
            unname(fit$econ_weights),
            unname(dat$w[expected]),
            tolerance = 1e-12
        )

        testthat::expect_equal(
            unname(fit$weights),
            unname(dat$w[expected]),
            tolerance = 1e-12
        )

        ###### 2. Construct systematic utilities ####

        X = stats::model.matrix(
            fit$terms,
            fit$econ_model,
            contrasts.arg = fit$contrasts
        )

        y = stats::model.response(
            fit$econ_model
        )

        eta = cbind(
            0,
            X %*% t(coef(fit))
        )

        colnames(eta) = fit$lev

        ###### 3. Compute raw observed-choice probabilities ####

        observed_prob = .econ_multinomial_probit_probs_general(
            eta = eta,
            Omega = fit$Omega,
            observed = as.integer(y)
        )

        ###### 4. Calculate independent weighted likelihood ####

        ll_reference = sum(
            dat$w[expected] *
                log(
                    pmax(
                        observed_prob,
                        .Machine$double.xmin
                    )
                )
        )

        ###### 5. Validate weighted log-likelihood ####

        testthat::expect_equal(
            as.numeric(logLik(fit)),
            ll_reference,
            tolerance = 1e-6
        )
    }
)


###### 49.4. Cluster Alignment After Sample Selection ####

testthat::test_that(
    "Estimated probit aligns clusters with the estimation sample",
    {

        z = .stage49_fixture()

        fit = z$fit
        dat = z$dat
        expected = z$expected

        testthat::expect_equal(
            unname(fit$cluster),
            unname(dat$cl[expected])
        )

        testthat::expect_equal(
            length(fit$cluster),
            nobs(fit)
        )

        testthat::expect_false(
            anyNA(fit$cluster)
        )

        testthat::expect_identical(
        colnames(fit$vcov),
        fit$parameter_names
    )    }
)


###### 49.5. Missing Weight Validation ####

testthat::test_that(
    "Estimated probit drops missing observation weights",
    {

        ###### 1. Generate data ####

        dat = .stage18_fixture()$dat

        dat$w = rep(1, nrow(dat))
        dat$w[5] = NA_real_

        ###### 2. Missing weights without subset ####

        testthat::expect_warning(
            {
                fit = econ_multinomial(
                    y ~ x,
                    data = dat,
                    weights = w,
                    link = "probit",
                    error_covariance = "estimated"
                )
            },
            "missingness in the weights"
        )

        expected = which(!is.na(dat$w))

        testthat::expect_equal(
            fit$used,
            expected
        )

        testthat::expect_equal(
            fit$omitted,
            5L
        )

        testthat::expect_equal(
            nobs(fit),
            length(expected)
        )

        testthat::expect_equal(
            unname(fit$econ_weights),
            unname(dat$w[expected])
        )

        ###### 3. Missing weights with subset ####

        dat$keep = rep(TRUE, nrow(dat))
        dat$keep[5] = FALSE

        testthat::expect_warning(
    {
        fit_subset = econ_multinomial(
            y ~ x,
            data = dat,
            weights = w,
            subset = keep,
            link = "probit",
            error_covariance = "estimated"
        )
    },
    "missingness in the weights"
)

        testthat::expect_equal(
            fit_subset$used,
            which(dat$keep)
        )

        testthat::expect_equal(
            nobs(fit_subset),
            sum(dat$keep)
        )

        testthat::expect_false(
            anyNA(fit_subset$econ_weights)
        )

    }
)

###### 49.6. Invalid Weight Validation ####

testthat::test_that(
    "Estimated probit rejects invalid observation weights",
    {

        dat = .stage18_fixture()$dat

        dat$w = rep(1, nrow(dat))

        for (invalid in c(-1, Inf, -Inf, NaN)) {

            dat$w[1] = invalid

            testthat::expect_error(
                econ_multinomial(
                    y ~ x,
                    data = dat,
                    weights = w,
                    link = "probit",
                    error_covariance = "estimated"
                ),
                "weights"
            )
        }
    }
)


###### 49.7. Out-of-Sample Missing Predictor Handling ####

testthat::test_that(
    "Estimated probit preserves missing values in out-of-sample predictions",
    {

        fit = .stage49_fixture()$fit

        newdata = data.frame(
            x = c(-1, NA_real_, 0, 1)
        )

        P = predict(
            fit,
            newdata = newdata,
            type = "probs"
        )

        testthat::expect_equal(
            dim(P),
            c(4L, 3L)
        )

        testthat::expect_true(
            all(is.na(P[2L, ]))
        )

        testthat::expect_true(
            all(is.finite(P[c(1L, 3L, 4L), ]))
        )

        testthat::expect_equal(
            rowSums(P[c(1L, 3L, 4L), ]),
            rep(1, 3L),
            tolerance = 1e-7
        )

        predicted_class = predict(
            fit,
            newdata = newdata,
            type = "class"
        )

        testthat::expect_true(
            is.na(predicted_class[2L])
        )
    }
)


###### 49.8. Weighted Score Contributions ####

testthat::test_that(
    "Estimated probit scores include observation weights",
    {

        z = .stage49_fixture()

        fit = z$fit

        X = stats::model.matrix(
            fit$terms,
            fit$econ_model,
            contrasts.arg = fit$contrasts
        )

        y = stats::model.response(
            fit$econ_model
        )

        selected = unique(as.integer(round(seq(
        1, nrow(X), length.out = min(36L, nrow(X))
    ))))

    X = X[selected, , drop = FALSE]
    y = y[selected]
    weights_selected = fit$econ_weights[selected]

    S_weighted = .econ_multinomial_probit_scores(
        fit = fit,
        X = X,
        y = y,
        weights = weights_selected
    )

    S_unweighted = .econ_multinomial_probit_scores(
            fit = fit,
            X = X,
            y = y,
            weights = rep(1, nrow(X))
        )

        reference = sweep(
            S_unweighted,
            1L,
            weights_selected,
            `*`
        )

        testthat::expect_equal(
            unname(S_weighted),
            unname(reference),
            tolerance = 1e-7
        )

        testthat::expect_identical(
            colnames(S_weighted),
            colnames(fit$vcov)
        )

        testthat::expect_true(
            all(is.finite(S_weighted))
        )
    }
)