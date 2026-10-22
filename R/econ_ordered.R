###### 1. Ordered Logit and Probit Estimation ####

#' Ordered Logit and Probit Regression
#'
#' @param formula Model formula.
#' @param data Data frame containing model variables.
#' @param link Either "logit" or "probit".
#' @param weights Optional observation weights.
#' @param subset Optional subset of observations.
#' @param na.action Missing-value handling function.
#' @param vcov Covariance estimator: "classical", "HC0", "HC1",
#'   "CR0", "CR1", "CR1p", "CR1S", or "stata".
#' @param cluster Optional clustering variable.
#' @return An object of class "econ_ordered" and "polr".
#' @export
econ_ordered = function(
    formula,
    data,
    link = c("logit", "probit"),
    weights = NULL,
    subset = NULL,
    na.action = stats::na.omit,
    vcov = c(
        "classical",
        "HC0", "HC1",
        "CR0", "CR1", "CR1p", "CR1S", "stata"),
    cluster = NULL
) {

    link = match.arg(link)
    vcov = match.arg(vcov)

    n = nrow(data)

    # Evaluate weights without prematurely forcing the argument
    weights_expr = substitute(weights)

    if (missing(weights) || identical(weights_expr, quote(NULL))) {

        evaluated_weights = rep(1, n)

    } else {

        evaluated_weights = eval(
            weights_expr,
            envir = data,
            enclos = parent.frame()
        )

        if (length(evaluated_weights) != n) {
            stop("weights must have one value per row of data.")
        }
    }

    # Evaluate subset in the original data
    subset_expr = substitute(subset)

    if (missing(subset) || identical(subset_expr, quote(NULL))) {

        evaluated_subset = rep(TRUE, n)

    } else {

        evaluated_subset = eval(
            subset_expr,
            envir = data,
            enclos = parent.frame()
        )

        if (length(evaluated_subset) != n) {
            stop("subset must have one value per row of data.")
        }

        if (!is.logical(evaluated_subset)) {
            stop("subset must be a logical vector.")
        }
    }

    # Evaluate cluster in the original data
    cluster_expr = substitute(cluster)

    if (missing(cluster) || identical(cluster_expr, quote(NULL))) {

        evaluated_cluster = NULL

    } else {

        evaluated_cluster = eval(
            cluster_expr,
            envir = data,
            enclos = parent.frame()
        )

        if (length(evaluated_cluster) != n) {
            stop("cluster must have one value per row of data.")
        }
    }

###### 1. Validate clustered covariance requirements ####

clustered_types = c(
    "CR0", "CR1", "CR1p", "CR1S", "stata"
)

if (vcov %in% clustered_types && is.null(evaluated_cluster)) {
    stop(
        "'", vcov, "' covariance estimation requires a cluster variable.",
        call. = FALSE
    )
}

    # Construct model frame
    model_data = data

    model_data$.econ_ordered_weights = evaluated_weights
    model_data$.econ_ordered_subset = evaluated_subset

    mf = stats::model.frame(
        formula = formula,
        data = model_data,
        weights = .econ_ordered_weights,
        subset = .econ_ordered_subset,
        na.action = na.action,
        drop.unused.levels = TRUE
    )

    model_weights = stats::model.weights(mf)

    if (anyNA(model_weights) ||
        any(!is.finite(model_weights)) ||
        any(model_weights < 0)) {

        stop("weights must be finite, nonnegative, and nonmissing.")
    }

    # Validate ordered response
    y = stats::model.response(mf)

    if (!is.ordered(y)) {
        stop("The dependent variable must be an ordered factor.")
    }

    if (nlevels(y) < 3L) {
        stop("The dependent variable must have at least three ordered levels.")
    }

    # Prepare estimation data
    estimation_data = mf

    estimation_data[["(.econ_ordered_weights)"]] = NULL

    estimation_data$.econ_prior_weights = model_weights

    # Estimate ordered model using MASS
    fit = MASS::polr(
        formula = formula,
        data = estimation_data,
        weights = .econ_prior_weights,
        method = if (link == "logit") "logistic" else "probit",
        Hess = TRUE,
        model = TRUE
    )

    # Store EconR information
    fit$na.action = attr(mf, "na.action")
    fit$link = link
    fit$vcov_type = vcov
    fit$weights = model_weights

    fit$used = as.integer(
        match(rownames(mf), rownames(data))
    )

    fit$omitted = setdiff(seq_len(n), fit$used)

    # 2. Align and validate cluster identifiers #

    if (is.null(evaluated_cluster)) {

        fit$cluster = NULL

    } else {

    fit$cluster = evaluated_cluster[fit$used]

    if (vcov %in% clustered_types) {

        if (anyNA(fit$cluster)) {
            stop(
                "'cluster' must not contain missing values.",
                call. = FALSE
            )
        }

        G = length(unique(fit$cluster))

        if (G < 2L) {
            stop(
                "'", vcov,
                "' covariance estimation requires at least two distinct clusters.",
                call. = FALSE
            )
        }
    }
}

    class(fit) = c("econ_ordered", class(fit))

    return(fit)
}


###### 2. Ordered Model Score Contributions ####

#' @export
estfun.econ_ordered = function(x, ...) {

    X = stats::model.matrix(x)

    X = X[
        ,
        colnames(X) != "(Intercept)",
        drop = FALSE
    ]

    beta = stats::coef(x)
    zeta = x$zeta

    y = as.integer(stats::model.response(x$model))
    w = x$weights

    eta = drop(X %*% beta)

    n = length(y)
    K = length(beta)
    J = length(zeta) + 1L

    if (x$link == "logit") {

        F = stats::plogis
        f = stats::dlogis

    } else {

        F = stats::pnorm
        f = stats::dnorm
    }

    scores = matrix(
        0,
        nrow = n,
        ncol = K + J - 1L
    )

    colnames(scores) = c(names(beta), names(zeta))

    for (i in seq_len(n)) {

        j = y[i]

        upper = if (j == J) Inf else zeta[j]
        lower = if (j == 1L) -Inf else zeta[j - 1L]

        a = upper - eta[i]
        b = lower - eta[i]

        p = F(a) - F(b)

        f_upper = if (j == J) 0 else f(a)
        f_lower = if (j == 1L) 0 else f(b)

        # Regression coefficient scores
        scores[i, seq_len(K)] =
            w[i] * X[i, ] * (f_lower - f_upper) / p

        # Threshold scores
        if (j < J) {
            scores[i, K + j] = w[i] * f_upper / p
        }

        if (j > 1L) {
            scores[i, K + j - 1L] =
                -w[i] * f_lower / p
        }
    }

    return(scores)
}


###### 3. Ordered Model Covariance Estimation ####

#' @export
vcov.econ_ordered = function(object, ...) {

    ###### 1. Classical covariance ####

    V_classical = NextMethod()

    if (object$vcov_type == "classical") {
        return(V_classical)
    }

    ###### 2. Individual score contributions ####

    S = sandwich::estfun(object)

    n = nrow(S)
    k = ncol(S)

    ###### 3. Heteroskedasticity-robust covariance ####

    if (object$vcov_type %in% c("HC0", "HC1")) {

    M = crossprod(S)

    if (object$vcov_type == "HC1") {

        if (n <= k) {
            stop(
                "HC1 requires more observations than estimated parameters.",
                call. = FALSE
            )
        }

        M = M * n / (n - k)

    } 

    ###### 4. Cluster-robust covariance ####

    } else if (object$vcov_type %in% c(
        "CR0", "CR1", "CR1p", "CR1S", "stata"
    )) {

        if (is.null(object$cluster)) {
            stop(
                "Clustered covariance requires a cluster variable.",
                call. = FALSE
            )
        }

        cluster = as.factor(object$cluster)

        if (length(cluster) != n) {
            stop(
                "Cluster variable is not aligned with the estimation sample.",
                call. = FALSE
            )
        }

        if (anyNA(cluster)) {
            stop(
                "Cluster variable contains missing values.",
                call. = FALSE
            )
        }

        G = nlevels(droplevels(cluster))

        if (G < 2L) {
            stop(
                "Clustered covariance requires at least two clusters.",
                call. = FALSE
            )
        }

        cluster_scores = rowsum(
            S,
            group = cluster,
            reorder = TRUE
        )

        M = crossprod(cluster_scores)

        ###### 5. Finite-sample corrections ####

        if (object$vcov_type == "CR1") {

            # Preserve existing EconR CR1 behavior

            if (n <= k) {
                stop("CR1 requires n > k.", call. = FALSE)
            }

            M = M * (n - 1) / (n - k)

        } else if (object$vcov_type == "CR1p") {

            if (G <= k) {
                stop(
                    "CR1p requires more clusters than parameters (G > k).",
                    call. = FALSE
                )
            }

            M = M * G / (G - k)

        } else if (object$vcov_type %in% c("CR1S", "stata")) {

            if (n <= k) {
                stop(
                    "CR1S and stata require n > k.",
                    call. = FALSE
                )
            }

            M = M * (
                G * (n - 1) /
                ((G - 1) * (n - k))
            )
        }

    } else {

        stop(
            "Unsupported covariance estimator: ",
            object$vcov_type,
            call. = FALSE
        )
    }

    ###### 6. Sandwich covariance ####

    V = V_classical %*% M %*% V_classical

    dimnames(V) = dimnames(V_classical)

    return(V)
}

###### 3a. McFadden Pseudo-R-Squared ####

.econ_ordered_pseudo_r2 = function(object) {

    ###### 1. Extract Outcome and Weights ####

    y = stats::model.response(object$model)
    w = object$weights

    ###### 2. Calculate Weighted Category Counts ####

    counts = as.numeric(
        rowsum(
            w,
            group = as.integer(y),
            reorder = TRUE
        )
    )

    counts = counts[counts > 0]

    if (length(counts) < 2L) {
        return(NA_real_)
    }

    ###### 3. Calculate Null Log-Likelihood ####

    probabilities = counts / sum(counts)

    null_loglik = sum(
        counts * log(probabilities)
    )

    ###### 4. Calculate McFadden Pseudo-R-Squared ####

    fitted_loglik = as.numeric(stats::logLik(object))

    if (!is.finite(fitted_loglik) ||
        !is.finite(null_loglik) ||
        null_loglik >= 0) {
        return(NA_real_)
    }

    1 - fitted_loglik / null_loglik
}

###### 4. Ordered Model Summary ####

#' @export
summary.econ_ordered = function(
    object,
    extra.significance = FALSE,
    exponentiate = FALSE,
    ...
) {

    ###### 1. Extract coefficients and covariance ####

    estimate = c(
        stats::coef(object),
        object$zeta
    )

    V = stats::vcov(object)

    std_error = sqrt(diag(V))
    std_error = std_error[names(estimate)]


        ###### 2. Determine inference distribution ####

    clustered_types = c(
        "CR0", "CR1", "CR1p", "CR1S", "stata"
    )

    is_clustered = object$vcov_type %in% clustered_types

    if (is_clustered) {

        G = length(unique(object$cluster))

        if (G < 2L) {
            stop(
                "Clustered inference requires at least two clusters.",
                call. = FALSE
            )
        }

        df = G - 1L

    } else {

        df = Inf
    }


    ###### 3. Compute coefficient inference ####

    statistic = estimate / std_error

    if (is_clustered) {

        p_value = 2 * stats::pt(
            abs(statistic),
            df = df,
            lower.tail = FALSE
        )

        statistic_name = "t value"
        p_value_name = "Pr(>|t|)"

    } else {

        p_value = 2 * stats::pnorm(
            abs(statistic),
            lower.tail = FALSE
        )

        statistic_name = "z value"
        p_value_name = "Pr(>|z|)"
    }


    ###### 4. Construct coefficient table ####

    coef_table = cbind(
        Estimate = estimate,
        "Std. Error" = std_error,
        statistic
    )

    colnames(coef_table)[3L] = statistic_name

    coef_table = cbind(
        coef_table,
        p_value
    )

    colnames(coef_table)[4L] = p_value_name

        ###### 5. Exponentiate coefficients ####

    if (!is.logical(exponentiate) ||
        length(exponentiate) != 1L ||
        is.na(exponentiate)) {

        stop("'exponentiate' must be TRUE or FALSE.", call. = FALSE)
    }

    if (exponentiate && object$link != "logit") {
        stop(
            "Exponentiation is only supported for ordered Logit models.",
            call. = FALSE
        )
    }

    if (exponentiate) {

        is_threshold = rownames(coef_table) %in% names(object$zeta)

        coef_table[!is_threshold, "Estimate"] =
            exp(coef_table[!is_threshold, "Estimate"])

        coef_table[is_threshold, "Estimate"] =
            exp(-coef_table[is_threshold, "Estimate"])
    }


    ###### 5a. Construct summary object ####

    ans = list(
        call = object$call,
        coefficients = coef_table,
        link = object$link,
        vcov_type = object$vcov_type,
        inference_distribution = if (is_clustered) "t" else "z",
        inference_df = df,
        residual.deviance = object$deviance,
        AIC = object$deviance + 2 * length(estimate),
        pseudo.r.squared = .econ_ordered_pseudo_r2(object),
        df.residual = object$df.residual,
        extra.significance = extra.significance,
        exponentiate = exponentiate
    )

    class(ans) = "summary.econ_ordered"

    ans
}


###### 5. Print Ordered Model Summary ####

#' @export
print.summary.econ_ordered = function(
    x,
    ...
) {

    ###### 1. Print call ####

    cat("\nCall:\n")

    print(x$call)


    ###### 2. Print model information ####

    cat(
        "\nOrdered response model: ",
        x$link,
        "\nCovariance type: ",
        x$vcov_type,
        "\n",
        sep = ""
    )


    ###### 3. Print coefficients ####

    ###### 3a. Print inference information ####

    if (identical(x$inference_distribution, "t")) {

        cat(
            "Inference: t distribution (df = ",
            x$inference_df,
            ")\n",
            sep = ""
        )

    } else {

        cat("Inference: standard normal distribution\n")
    }

    cat("\nCoefficients:\n")

        if (isTRUE(x$exponentiate)) {
        cat(
            "Odds ratios (standard errors and test statistics ",
            "on original coefficient scale):\n",
            sep = ""
        )
    }

    coef_table = x$coefficients

    ##### Significance symbols ####

    p_value_column = if (identical(x$inference_distribution, "t")) {
    "Pr(>|t|)"
} else {
    "Pr(>|z|)"
}

p_values = coef_table[, p_value_column]
    
    if (isTRUE(x$extra.significance)) {

        stars = ifelse(
            is.na(p_values), "",
            ifelse(
                p_values < 0.001, "+",
                ifelse(
                    p_values < 0.01, "***",
                    ifelse(
                        p_values < 0.05, "**",
                        ifelse(p_values < 0.10, "*", "")
                    )
                )
            )
        )

    } else {

        stars = ifelse(
            is.na(p_values), "",
            ifelse(
                p_values < 0.01, "***",
                ifelse(
                    p_values < 0.05, "**",
                    ifelse(p_values < 0.10, "*", "")
                )
            )
        )
    }


    ###### 4. Format coefficient table ####

    formatted = matrix(
        "",
        nrow = nrow(coef_table),
        ncol = ncol(coef_table),
        dimnames = dimnames(coef_table)
    )

    for (j in seq_len(ncol(coef_table))) {

        formatted[, j] = sprintf(
            "%.5f",
            coef_table[, j]
        )
    }


    ###### 5. Add significance symbols ####

    formatted = cbind(
        formatted,
        " " = stars
    )

    print(
        formatted,
        quote = FALSE,
        right = TRUE
    )


    ###### 6. Print significance legend ####

    if (isTRUE(x$extra.significance)) {

        cat(
            "\nSignificance: * p < 0.1, ** p < 0.05, *** p < 0.01, + p < 0.001\n"
        )

    } else {

        cat(
            "\nSignificance: * p < 0.1, ** p < 0.05, *** p < 0.01\n"
        )
    }


    ###### 7. Print model fit ####

    cat(
        "\nResidual deviance: ",
        sprintf("%.5f", x$residual.deviance),
        "\nAIC: ",
        sprintf("%.5f", x$AIC),
        "\n",
        sep = ""
    )

    invisible(x)

    ###### 8. Print McFadden Pseudo-R-Squared ####

    cat(
        "McFadden pseudo-R-squared: ",
        sprintf("%.5f", x$pseudo.r.squared),
        "\n",
        sep = ""
    )

}


###### 6. Global Variables ####

utils::globalVariables(c(
    ".econ_ordered_weights",
    ".econ_ordered_subset"
))


###### 7. Ordered Model Tidy Method ####

#' @export
tidy.econ_ordered = function(
    x,
    conf.int = FALSE,
    conf.level = 0.95,
    exponentiate = FALSE,
    ...
) {

    ###### 1. Extract model summary ####

    s = summary(x)

    tab = s$coefficients

    ###### Validate exponentiation ####

    if (!is.logical(exponentiate) ||
        length(exponentiate) != 1L ||
        is.na(exponentiate)) {

        stop("'exponentiate' must be TRUE or FALSE.", call. = FALSE)
    }

    if (exponentiate && x$link != "logit") {
        stop(
            "Exponentiation is only supported for ordered Logit models.",
            call. = FALSE
        )
    }


    ###### 2. Determine inference distribution ####

    is_clustered = identical(
        s$inference_distribution,
        "t"
    )

    statistic_column = if (is_clustered) {
        "t value"
    } else {
        "z value"
    }

    p_value_column = if (is_clustered) {
        "Pr(>|t|)"
    } else {
        "Pr(>|z|)"
    }


    ###### 3. Construct tidy results ####

    result = data.frame(
        term = rownames(tab),
        estimate = unname(tab[, "Estimate"]),
        std.error = unname(tab[, "Std. Error"]),
        statistic = unname(tab[, statistic_column]),
        p.value = unname(tab[, p_value_column]),
        coef.type = ifelse(
            rownames(tab) %in% names(x$zeta),
            "scale",
            "coefficient"
        ),
        row.names = NULL
    )


    ###### 4. Calculate confidence intervals ####

    if (isTRUE(conf.int)) {

        if (!is.numeric(conf.level) ||
            length(conf.level) != 1L ||
            is.na(conf.level) ||
            !is.finite(conf.level) ||
            conf.level <= 0 ||
            conf.level >= 1) {

            stop(
                "'conf.level' must be a number strictly between 0 and 1.",
                call. = FALSE
            )
        }

        alpha = 1 - conf.level

        if (is_clustered) {

            critical = stats::qt(
                1 - alpha / 2,
                df = s$inference_df
            )

        } else {

            critical = stats::qnorm(
                1 - alpha / 2
            )
        }

        result$conf.low = result$estimate -
            critical * result$std.error

        result$conf.high = result$estimate +
            critical * result$std.error
    }


        ###### 5. Exponentiate estimates and confidence intervals ####

    if (exponentiate) {

        is_threshold = result$coef.type == "scale"

        if (isTRUE(conf.int)) {

            lower = result$conf.low
            upper = result$conf.high

            result$conf.low = exp(lower)
            result$conf.high = exp(upper)

            result$conf.low[is_threshold] = exp(-upper[is_threshold])
            result$conf.high[is_threshold] = exp(-lower[is_threshold])
        }

        result$estimate[!is_threshold] =
            exp(result$estimate[!is_threshold])

        result$estimate[is_threshold] =
            exp(-result$estimate[is_threshold])
    }

    ###### 6. Return tidy results ####

    tibble::as_tibble(result)
}


###### 8. Ordered Model Confidence Intervals ####

#' @export
confint.econ_ordered = function(
    object,
    parm,
    level = 0.95,
    ...
) {

    ###### 1. Validate confidence level ####

    if (!is.numeric(level) ||
        length(level) != 1L ||
        is.na(level) ||
        !is.finite(level) ||
        level <= 0 ||
        level >= 1) {

        stop(
            "'level' must be a single number strictly between 0 and 1.",
            call. = FALSE
        )
    }


    ###### 2. Extract all model parameters ####

    b = c(
        stats::coef(object),
        object$zeta
    )

    se = sqrt(diag(stats::vcov(object)))

    se = se[names(b)]


    ###### 3. Select parameters ####

    if (missing(parm)) {

        parm = seq_along(b)

    } else if (is.character(parm)) {

        parm = match(parm, names(b))

        if (anyNA(parm)) {
            stop(
                "Unknown parameter name in 'parm'.",
                call. = FALSE
            )
        }

    } else if (is.numeric(parm)) {

        if (anyNA(parm) ||
            any(!is.finite(parm)) ||
            any(parm != floor(parm)) ||
            any(parm < 1 | parm > length(b))) {

            stop(
                "Invalid parameter selection in 'parm'.",
                call. = FALSE
            )
        }

    } else {

        stop(
            "'parm' must contain parameter names or indices.",
            call. = FALSE
        )
    }


    ###### 4. Determine inference distribution ####

    s = summary(object)

    if (identical(s$inference_distribution, "t")) {

        critical = stats::qt(
            (1 + level) / 2,
            df = s$inference_df
        )

    } else {

        critical = stats::qnorm(
            (1 + level) / 2
        )
    }


    ###### 5. Construct confidence intervals ####

    ci = cbind(
        b[parm] - critical * se[parm],
        b[parm] + critical * se[parm]
    )

    alpha = (1 - level) / 2

    colnames(ci) = paste0(
        format(
            100 * c(alpha, 1 - alpha),
            trim = TRUE
        ),
        " %"
    )

    ci
}



#### 9. Ordered Model ANOVA ####

#' @export
anova.econ_ordered = function(object, ..., test = "Chisq") {

    ###### 1. Collect and validate models ####

    models = c(list(object), list(...))

    if (length(models) < 2L) {
        stop(
            "Ordered ANOVA currently requires at least two models.",
            call. = FALSE
        )
    }

    if (!all(vapply(
        models,
        inherits,
        logical(1),
        what = "econ_ordered"
    ))) {
        stop(
            "All models must be econ_ordered objects.",
            call. = FALSE
        )
    }

    if (!identical(test, "Chisq")) {
        stop(
            "Only likelihood-ratio tests ('Chisq') are supported.",
            call. = FALSE
        )
    }

    ###### 2. Verify comparable estimation data ####

    reference = models[[1L]]

    for (i in seq_along(models)[-1L]) {

        current = models[[i]]

        if (!identical(reference$link, current$link)) {
            stop("Models must use the same link function.", call. = FALSE)
        }

        if (!identical(reference$used, current$used) ||
            !identical(
                rownames(reference$model),
                rownames(current$model)
            )) {
            stop(
                "Models must use identical estimation samples.",
                call. = FALSE
            )
        }

        y_ref = stats::model.response(reference$model)
        y_cur = stats::model.response(current$model)

        if (!identical(y_ref, y_cur)) {
            stop(
                "Models must use identical dependent-variable values.",
                call. = FALSE
            )
        }

                if (!isTRUE(all.equal(
            reference$weights,
            current$weights,
            check.attributes = FALSE
        ))) {
            stop(
                "Models must use identical observation weights.",
                call. = FALSE
            )
        }

        ###### Verify shared regressor values ####

        X_ref = stats::model.matrix(reference)
        X_cur = stats::model.matrix(current)

        shared = intersect(
            colnames(X_ref),
            colnames(X_cur)
        )

        if (length(shared) > 0L) {
            if (!isTRUE(all.equal(
                unname(X_ref[, shared, drop = FALSE]),
                unname(X_cur[, shared, drop = FALSE]),
                tolerance = 1e-10
            ))) {
                stop(
                    "Models contain different values for shared regressors.",
                    call. = FALSE
                )
            }
        }
    }

    ###### 3. Extract likelihoods and model dimensions ####

    ll = vapply(
        models,
        function(x) as.numeric(stats::logLik(x)),
        numeric(1)
    )

    df = vapply(
        models,
        function(x) as.numeric(attr(stats::logLik(x), "df")),
        numeric(1)
    )

    ###### 4. Verify nested model specifications ####

    X = lapply(models, function(x) {
        stats::model.matrix(x)
    })

    for (i in seq_len(length(models) - 1L)) {

        for (j in seq.int(i + 1L, length(models))) {

            small = if (df[i] <= df[j]) i else j
            large = if (df[i] <= df[j]) j else i

            if (df[small] == df[large]) {
                stop(
                    "Models must have different parameter counts.",
                    call. = FALSE
                )
            }

            rank_large = qr(X[[large]])$rank
            rank_combined = qr(
                cbind(X[[large]], X[[small]])
            )$rank

            if (rank_combined > rank_large) {
                stop(
                    "Model specifications must be nested.",
                    call. = FALSE
                )
            }
        }
    }

    ###### 5. Compute likelihood-ratio comparisons ####

    delta_df = c(NA_real_, diff(df))
    delta_ll = c(NA_real_, diff(ll))

    chisq = 2 * delta_ll
    chisq[1L] = NA_real_

    # Negative differences in model complexity represent reversed
    # comparisons, not valid forward likelihood-ratio tests.
    if (any(delta_df[-1L] <= 0)) {
        stop(
            "Models must be supplied in increasing complexity.",
            call. = FALSE
        )
    }

    if (any(chisq[-1L] < -1e-7)) {
        stop(
            "The larger model has a lower maximized log-likelihood.",
            call. = FALSE
        )
    }

    chisq[-1L] = pmax(chisq[-1L], 0)

    p_value = stats::pchisq(
        chisq,
        df = delta_df,
        lower.tail = FALSE
    )

    ###### 6. Construct ANOVA table ####

    result = data.frame(
        Model = vapply(
            models,
            function(x) paste(deparse(stats::formula(x)), collapse = ""),
            character(1)
        ),
        Df = df,
        logLik = ll,
        `Df Diff` = delta_df,
        Chisq = chisq,
        `Pr(>Chisq)` = p_value,
        check.names = FALSE
    )

    rownames(result) = seq_along(models)

    structure(
        result,
        heading = "Likelihood ratio tests of ordered regression models",
        class = c("anova", "data.frame")
    )
}



#### 10. Ordered Model Augmentation ####

#' @export
augment.econ_ordered = function(
    x,
    data = NULL,
    newdata = NULL,
    se_fit = FALSE,
    ...
) {

    ###### 1. Validate arguments ####

    if (!is.logical(se_fit) ||
        length(se_fit) != 1L ||
        is.na(se_fit)) {
        stop(
            "'se_fit' must be TRUE or FALSE.",
            call. = FALSE
        )
    }


    if (!is.null(data) && !is.null(newdata)) {
        stop(
            "Supply either 'data' or 'newdata', not both.",
            call. = FALSE
        )
    }

    ###### 2. Determine prediction data ####

    out_of_sample = !is.null(newdata)

    if (out_of_sample) {

        prediction_data = as.data.frame(newdata)
        result = prediction_data

    } else {

        prediction_data = as.data.frame(x$model)

        if (is.null(data)) {
            result = as.data.frame(x$model)
        } else {
            result = as.data.frame(data)
        }

        if (nrow(result) != length(x$used)) {
            stop(
                "'data' must contain the estimation observations.",
                call. = FALSE
            )
        }

        # Remove internal estimation variables
        internal = c(
            ".econ_prior_weights",
            ".econ_ordered_weights",
            ".econ_ordered_subset"
        )

        result = result[
            ,
            !names(result) %in% internal,
            drop = FALSE
        ]
    }

    ###### 3. Predict category probabilities ####

    probabilities = stats::predict(
        x,
        newdata = prediction_data,
        type = "probs"
    )

    probabilities = as.matrix(probabilities)

    categories = colnames(probabilities)

    if (is.null(categories)) {
        stop(
            "Predicted probabilities must have category names.",
            call. = FALSE
        )
    }

    ###### 4. Predict most likely category ####

    fitted_class = stats::predict(
        x,
        newdata = prediction_data,
        type = "class"
    )

    result$.fitted = fitted_class

    ###### 5. Add category probabilities ####

    probability_names = paste0(
        ".fitted_",
        make.names(categories, unique = TRUE)
    )

    for (j in seq_along(probability_names)) {
        result[[probability_names[j]]] =
            unname(probabilities[, j])
    }
    
    ###### 5a. Delta-method prediction standard errors ####

    if (isTRUE(se_fit)) {

        terms_object = stats::delete.response(stats::terms(x))
        mf_new = stats::model.frame(
            terms_object,
            data = prediction_data,
            xlev = x$xlevels,
            na.action = stats::na.pass
        )

        X = stats::model.matrix(
            terms_object,
            data = mf_new,
            contrasts.arg = x$contrasts
        )

        X = X[, names(stats::coef(x)), drop = FALSE]

        beta = stats::coef(x)
        zeta = x$zeta

        eta = drop(X %*% beta)

        K = length(beta)
        J = length(zeta) + 1L

        V = stats::vcov(x)
        parameter_names = c(names(beta), names(zeta))
        V = V[parameter_names, parameter_names, drop = FALSE]

        density = if (x$link == "logit") {
            stats::dlogis
        } else {
            stats::dnorm
        }

        se_prob = matrix(
            NA_real_,
            nrow = nrow(X),
            ncol = J,
            dimnames = list(NULL, categories)
        )

        for (i in seq_len(nrow(X))) {

            for (j in seq_len(J)) {

                f_upper = if (j == J) {
                    0
                } else {
                    density(zeta[j] - eta[i])
                }

                f_lower = if (j == 1L) {
                    0
                } else {
                    density(zeta[j - 1L] - eta[i])
                }

                gradient = numeric(K + J - 1L)

                # Derivative with respect to slope coefficients
                gradient[seq_len(K)] =
                    (f_lower - f_upper) * X[i, ]

                # Derivatives with respect to thresholds
                if (j < J) {
                    gradient[K + j] = f_upper
                }

                if (j > 1L) {
                    gradient[K + j - 1L] = -f_lower
                }

                variance = drop(
                    crossprod(gradient, V %*% gradient)
                )

                se_prob[i, j] = sqrt(pmax(variance, 0))
            }
        }

        se_names = paste0(
            ".se.fit_",
            make.names(categories, unique = TRUE)
        )

        for (j in seq_along(se_names)) {
            result[[se_names[j]]] = se_prob[, j]
        }
    }


    ###### 6. Calculate in-sample residuals ####

    if (!out_of_sample) {

        observed = stats::model.response(x$model)

        # Expected ordinal category index
        expected_index = drop(
            probabilities %*% seq_along(categories)
        )

        # Ordinal-score residual:
        # observed category index minus expected index
        result$.resid = as.integer(observed) -
            expected_index
    }

    ###### 7. Return augmented data ####

    tibble::as_tibble(result)
}


#### 11. Ordered Model Residuals ####

#' @export
residuals.econ_ordered = function(
    object,
    type = c("response", "pearson"),
    ...
) {

    type = match.arg(type)

    y = as.integer(stats::model.response(object$model))

    probabilities = as.matrix(
        stats::predict(
            object,
            newdata = object$model,
            type = "probs"
        )
    )

    categories = seq_len(ncol(probabilities))

    expected = drop(probabilities %*% categories)

    response_residuals = y - expected

    if (type == "response") {
    return(stats::naresid(object$na.action, response_residuals))
}

    variance = drop(
        probabilities %*% (categories^2)
    ) - expected^2

    pearson_residuals = response_residuals /
    sqrt(pmax(variance, .Machine$double.eps))

stats::naresid(object$na.action, pearson_residuals)
}


#### 12. Ordered Model Glance Method ####

#' @export
glance.econ_ordered = function(x, ...) {

    ###### 1. Obtain Standard Ordered Model Statistics ####

    result = utils::getS3method(
        "glance",
        "polr"
    )(x, ...)

    ###### 2. Add McFadden Pseudo-R-Squared ####

    result$pseudo.r.squared =
        .econ_ordered_pseudo_r2(x)

    ###### 3. Return Model Statistics ####

    result
}