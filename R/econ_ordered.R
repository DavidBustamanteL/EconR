###### 1. Ordered Logit and Probit Estimation ####

#' Ordered Logit and Probit Regression
#'
#' @param formula Model formula.
#' @param data Data frame containing model variables.
#' @param link Either "logit" or "probit".
#' @param weights Optional observation weights.
#' @param subset Optional subset of observations.
#' @param na.action Missing-value handling function.
#' @param vcov Covariance estimator: "classical", "HC0", "HC1", "CR0", or "CR1".
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
    vcov = c("classical", "HC0", "HC1", "CR0", "CR1"),
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

    if (vcov %in% c("CR0", "CR1") && is.null(evaluated_cluster)) {
        stop("CR0 and CR1 require a cluster variable.")
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
    fit$link = link
    fit$vcov_type = vcov
    fit$weights = model_weights

    fit$used = as.integer(
        match(rownames(mf), rownames(data))
    )

    fit$omitted = setdiff(seq_len(n), fit$used)

    # Align cluster identifiers with estimation observations
    if (is.null(evaluated_cluster)) {

        fit$cluster = NULL

    } else {

        fit$cluster = evaluated_cluster[fit$used]

        if (anyNA(fit$cluster)) {
            stop("cluster contains missing values in the estimation sample.")
        }

        if (length(unique(fit$cluster)) < 2L) {
            stop("Clustered covariance requires at least two clusters.")
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

    # Classical covariance in original threshold parameterization
    V_classical = NextMethod()

    if (object$vcov_type == "classical") {
        return(V_classical)
    }

    S = sandwich::estfun(object)

    if (object$vcov_type %in% c("HC0", "HC1")) {

        M = crossprod(S)

        if (object$vcov_type == "HC1") {

            n = nrow(S)
            k = ncol(S)

            if (n <= k) {
                stop(
                    "HC1 requires more observations than estimated parameters."
                )
            }

            M = M * n / (n - k)
        }

    } else if (object$vcov_type %in% c("CR0", "CR1")) {

        if (is.null(object$cluster)) {
            stop("CR0 and CR1 require a cluster variable.")
        }

        cluster = as.factor(object$cluster)

        if (length(cluster) != nrow(S)) {
            stop(
                "Cluster variable is not aligned with the estimation sample."
            )
        }

        if (anyNA(cluster)) {
            stop("Cluster variable contains missing values.")
        }

        cluster_scores = rowsum(
            S,
            group = cluster,
            reorder = TRUE
        )

        M = crossprod(cluster_scores)

        if (object$vcov_type == "CR1") {

            n = nrow(S)
            k = ncol(S)

            if (n <= k) {
                stop("CR1 requires more observations than estimated parameters.")
            }

            M = M * (n - 1) / (n - k)
        }

    } else {

        stop(
            "Unsupported covariance estimator: ",
            object$vcov_type
        )
    }

    V = V_classical %*% M %*% V_classical

    dimnames(V) = dimnames(V_classical)

    return(V)
}


###### 4. Ordered Model Summary ####

#' @export
summary.econ_ordered = function(
    object,
    extra.significance = FALSE,
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


    ###### 2. Compute coefficient inference ####

    z_value = estimate / std_error

    p_value = 2 * stats::pnorm(
        abs(z_value),
        lower.tail = FALSE
    )


    ###### 3. Construct coefficient table ####

    coef_table = cbind(
        Estimate = estimate,
        "Std. Error" = std_error,
        "z value" = z_value,
        "Pr(>|z|)" = p_value
    )


    ###### 4. Construct summary object ####

    ans = list(
        call = object$call,
        coefficients = coef_table,
        link = object$link,
        vcov_type = object$vcov_type,
        residual.deviance = object$deviance,
        AIC = object$deviance + 2 * length(estimate),
        df.residual = object$df.residual,
        extra.significance = extra.significance
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

    cat("\nCoefficients:\n")

    coef_table = x$coefficients

    ##### Significance symbols ####

    p_values = coef_table[, "Pr(>|z|)"]
    
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
}


###### 6. Global Variables ####

utils::globalVariables(c(
    ".econ_ordered_weights",
    ".econ_ordered_subset"
))


###### 7. Ordered Model Tidy Method ####

#' @export
tidy.econ_ordered = function(x, conf.int = FALSE,
                             conf.level = 0.95, ...) {

    s = summary(x)

    tab = s$coefficients

    result = data.frame(
        term = rownames(tab),
        estimate = unname(tab[, "Estimate"]),
        std.error = unname(tab[, "Std. Error"]),
        statistic = unname(tab[, "z value"]),
        p.value = unname(tab[, "Pr(>|z|)"]),
        coef.type = ifelse(
            rownames(tab) %in% names(x$zeta),
            "scale",
            "coefficient"
        ),
        row.names = NULL
    )

    if (isTRUE(conf.int)) {
        critical = stats::qnorm(
            1 - (1 - conf.level) / 2
        )

        result$conf.low = result$estimate -
            critical * result$std.error

        result$conf.high = result$estimate +
            critical * result$std.error
    }

    tibble::as_tibble(result)
}