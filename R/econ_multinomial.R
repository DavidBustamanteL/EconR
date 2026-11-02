###### 1. Multinomial Logit Estimation ####

#' Multinomial Logit Regression
#'
#' Estimates a baseline-category multinomial logit model using nnet::multinom.
#' The first response level is the reference category. Covariance estimates
#' use observation-level score contributions and the inverse Hessian.
#'
#' @param formula Model formula.
#' @param data Data frame containing model variables.
#' @param weights Optional nonnegative observation weights.
#' @param subset Optional logical or integer row subset.
#' @param na.action Missing-value handling function.
#' @param vcov One of "classical", "HC0", "HC1", "CR0", "CR1",
#'   "CR1p", "CR1S", or "stata".
#' @param cluster Optional clustering variable, required for CR estimators.
#' @param error_covariance For probit, either "fixed" (default) or "estimated".
#' @param ... Additional arguments passed to nnet::multinom().
#' @return An object inheriting from "econ_multinomial" and "multinom".
#' @export
econ_multinomial = function(
    formula, data, weights = NULL, subset = NULL,
    na.action = stats::na.omit,
    link = c("logit", "probit"),
    vcov = c("classical", "HC0", "HC1", "CR0", "CR1", "CR1p", "CR1S", "stata"),
    cluster = NULL,
    error_covariance = c("fixed", "estimated"),
    ...
) {
    link = match.arg(link)
    vcov = match.arg(vcov)
    error_covariance = match.arg(error_covariance)

    if (link != "probit" && error_covariance == "estimated") {
        stop(
            "'error_covariance = \"estimated\"' is only supported ",
            "when link = \"probit\".",
            call. = FALSE
        )
    }

    if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
    n_original = nrow(data)
    enclos = parent.frame()
    evaluate_arg = function(expr, label, default) {
        if (identical(expr, quote(NULL))) return(default)
        value = eval(expr, envir = data, enclos = enclos)
        if (length(value) != n_original) {
            stop(label, " must have one value per row of data.", call. = FALSE)
        }
        value
    }
    w = evaluate_arg(substitute(weights), "weights", rep(1, n_original))
    sub = if (missing(subset) || identical(substitute(subset), quote(NULL))) {
        rep(TRUE, n_original)
    } else {
        eval(substitute(subset), envir = data, enclos = enclos)
    }
    cl = evaluate_arg(substitute(cluster), "cluster", NULL)
    clustered = vcov %in% c("CR0", "CR1", "CR1p", "CR1S", "stata")
    if (clustered && is.null(cl)) stop("Clustered covariance requires 'cluster'.", call. = FALSE)
    if (anyNA(w) || any(!is.finite(w)) || any(w < 0)) {
        stop("weights must be finite, nonnegative, and nonmissing.", call. = FALSE)
    }
    model_data = data
    model_data$.econ_multi_row_id = seq_len(n_original)
    model_data$.econ_multi_weights = w
    model_data$.econ_multi_subset = sub
    mf = stats::model.frame(
        formula, data = model_data,
        weights = .econ_multi_weights,
        subset = .econ_multi_subset,
        na.action = na.action, drop.unused.levels = TRUE
    )
    y = stats::model.response(mf)
    if (!is.factor(y) || is.ordered(y) || nlevels(y) < 3L) {
        stop("Response must be an unordered factor with at least three observed levels.", call. = FALSE)
    }
    used = match(rownames(mf), rownames(data))
    if (anyNA(used)) stop("Could not align estimation rows to input data.", call. = FALSE)
    cl_used = if (is.null(cl)) NULL else cl[used]
    if (clustered && (anyNA(cl_used) || length(unique(cl_used)) < 2L)) {
        stop("Cluster identifiers must be nonmissing and contain at least two clusters.", call. = FALSE)
    }
  
    
    ###### Probit estimation branch ####

    if (link == "probit") {

        if (length(list(...)) > 0L) {
            stop(
                "Additional '...' arguments are not supported for probit.",
                call. = FALSE
            )
        }

        terms_object = stats::terms(mf)

        X = stats::model.matrix(
            terms_object,
            data = mf
        )

        fit = .econ_multinomial_probit_fit_general(
          X = X,
          y = y,
          weights = stats::model.weights(mf),
          error_covariance = error_covariance
        )
      
        fit$error_covariance = error_covariance
      
        fit$call = match.call()
        fit$link = link
        fit$vcov_type = vcov
        fit$cluster = cl_used
        fit$used = used
        fit$omitted = setdiff(seq_len(n_original), used)
        fit$econ_weights = stats::model.weights(mf)
        fit$econ_model = mf
        fit$econ_formula = formula

        fit$terms = terms_object
        fit$contrasts = attr(X, "contrasts")
        fit$xlevels = stats::.getXlevels(terms_object, mf)
        fit$lev = levels(y)

        class(fit) = "econ_multinomial"

        return(fit)
    }

  
    fit_data = mf
    fit_data[["(.econ_multi_weights)"]] = NULL
    fit_data$.econ_multi_fit_weights = stats::model.weights(mf)
    fit = nnet::multinom(
        formula, data = fit_data, weights = .econ_multi_fit_weights,
        Hess = TRUE, model = TRUE, trace = FALSE, ...
    )
    fit$call = match.call()
    fit$link = link
    fit$vcov_type = vcov
    fit$cluster = cl_used
    fit$used = used
    fit$omitted = setdiff(seq_len(n_original), used)
    fit$econ_weights = stats::model.weights(mf)
    fit$econ_model = mf
    fit$econ_formula = formula
    class(fit) = c("econ_multinomial", class(fit))
    fit
}

###### 2. Observation-Level Score Contributions ####


#' @export
estfun.econ_multinomial = function(x, ...) {

    X = stats::model.matrix(
        x$terms,
        x$econ_model,
        contrasts.arg = x$contrasts
    )

    y = stats::model.response(x$econ_model)

    if (identical(x$link, "probit")) {
        return(.econ_multinomial_probit_scores(
            fit = x,
            X = X,
            y = y
        ))
    }

    p = stats::fitted(x)

    if (is.null(dim(p))) {
        stop(
            "Expected a multinomial probability matrix.",
            call. = FALSE
        )
    }

    lev = x$lev

    if (is.null(lev)) {
        lev = levels(y)
    }

    nonbase = lev[-1L]
    w = x$econ_weights

    parts = lapply(nonbase, function(category) {
        sweep(
            X,
            1L,
            w * (as.numeric(y == category) - p[, category]),
            `*`
        )
    })

    S = do.call(cbind, parts)

    colnames(S) = unlist(
        lapply(nonbase, function(category) {
            paste0(category, ":", colnames(X))
        }),
        use.names = FALSE
    )

    S
}



###### 3. Covariance Estimation ####

#' @export
vcov.econ_multinomial = function(object, ...) {

    ###### 1. Multinomial probit covariance ####

    if (identical(object$link, "probit")) {

        X = stats::model.matrix(
            object$terms,
            object$econ_model,
            contrasts.arg = object$contrasts
        )

        y = stats::model.response(object$econ_model)

        return(.econ_multinomial_probit_vcov(
            fit = object,
            X = X,
            y = y,
            vcov = object$vcov_type,
            cluster = object$cluster
        ))
    }

    ###### 2. Multinomial logit covariance ####

    V = NextMethod()

    if (object$vcov_type == "classical") return(V)

    S = estfun.econ_multinomial(object)

    if (!identical(colnames(S), colnames(V))) {
        if (!setequal(colnames(S), colnames(V))) {
            stop("Score and Hessian parameter names differ.", call. = FALSE)
        }
        S = S[, colnames(V), drop = FALSE]
    }

    n = nrow(S)
    k = ncol(S)

    if (object$vcov_type %in% c("HC0", "HC1")) {

        M = crossprod(S)

        if (object$vcov_type == "HC1") {
            if (n <= k) stop("HC1 requires n > k.", call. = FALSE)
            M = M * n / (n - k)
        }

    } else {

        cluster = object$cluster

        if (is.null(cluster) || length(cluster) != n || anyNA(cluster)) {
            stop("Invalid or misaligned cluster variable.", call. = FALSE)
        }

        G = length(unique(cluster))

        if (G < 2L) stop("At least two clusters are required.", call. = FALSE)

        U = rowsum(S, group = cluster, reorder = TRUE)
        M = crossprod(U)

        if (object$vcov_type == "CR1") {

            if (n <= k) stop("CR1 requires n > k.", call. = FALSE)
            M = M * (n - 1) / (n - k)

        } else if (object$vcov_type == "CR1p") {

            if (G <= k) stop("CR1p requires G > k.", call. = FALSE)
            M = M * G / (G - k)

        } else if (object$vcov_type %in% c("CR1S", "stata")) {

            if (n <= k) stop("CR1S/stata requires n > k.", call. = FALSE)

            M = M * G * (n - 1) / ((G - 1) * (n - k))
        }
    }

    result = V %*% M %*% V

    dimnames(result) = dimnames(V)

    (result + t(result)) / 2
}


###### 3.1. Utility-Coefficient Covariance Block ####

.econ_multinomial_beta_vcov = function(object) {
    V = stats::vcov(object)
    C = stats::coef(object)
    if (is.null(dim(C))) C = matrix(C, nrow = 1L)
    beta_names = unlist(lapply(rownames(C), function(category) {
        paste0(category, ":", colnames(C))
    }), use.names = FALSE)
    if (is.null(rownames(C)) || !all(beta_names %in% colnames(V))) {
        stop("Utility coefficient and covariance names differ.", call. = FALSE)
    }
    V[beta_names, beta_names, drop = FALSE]
}


###### 4. Model Summary and Inference ####

#' @export
summary.econ_multinomial = function(
    object,
    extra.significance = FALSE,
    ...
) {

      if (!is.logical(extra.significance) ||
        length(extra.significance) != 1L ||
        is.na(extra.significance)) {

        stop(
            "'extra.significance' must be TRUE or FALSE.",
            call. = FALSE
        )
    }

    ###### 1. Extract coefficient inference ####

    C = stats::coef(object)

    if (is.null(dim(C))) {
        C = matrix(C, nrow = 1L)
    }

    V = .econ_multinomial_beta_vcov(object)

    se = matrix(
        sqrt(diag(V)),
        nrow = nrow(C),
        byrow = TRUE,
        dimnames = dimnames(C)
    )

    stat = C / se

    clustered = object$vcov_type %in% c(
        "CR0", "CR1", "CR1p", "CR1S", "stata"
    )

    df = if (clustered) {
        length(unique(object$cluster)) - 1L
    } else {
        Inf
    }

    p = if (is.finite(df)) {
        2 * stats::pt(-abs(stat), df = df)
    } else {
        2 * stats::pnorm(-abs(stat))
    }

    ###### 2. Joint Wald test ####

    beta = as.vector(t(C))
    names(beta) = colnames(V)

    tested = !grepl(":(Intercept)$", names(beta))

    b = beta[tested]

    V_test = V[
        names(b),
        names(b),
        drop = FALSE
    ]

    q = length(b)

    if (q > 0L) {

        W = as.numeric(
            crossprod(b, solve(V_test, b))
        )

        if (clustered) {

            wald_statistic = W / q
            wald_distribution = "F"
            wald_df1 = q
            wald_df2 = df

            wald_p_value = stats::pf(
                wald_statistic,
                df1 = wald_df1,
                df2 = wald_df2,
                lower.tail = FALSE
            )

        } else {

            wald_statistic = W
            wald_distribution = "Chisq"
            wald_df1 = q
            wald_df2 = NA_real_

            wald_p_value = stats::pchisq(
                W,
                df = q,
                lower.tail = FALSE
            )
        }

    } else {

        wald_statistic = NA_real_
        wald_distribution = if (clustered) "F" else "Chisq"
        wald_df1 = 0L
        wald_df2 = if (clustered) df else NA_real_
        wald_p_value = NA_real_
    }

    ###### 3. Construct summary object ####

    out = list(
        link = object$link,
        call = object$call,
        coefficients = C,
        standard.errors = se,
        statistic = stat,
        p.value = p,
        vcov = V,
        vcov_type = object$vcov_type,
        df = df,
        logLik = stats::logLik(object),
        nobs = stats::nobs(object),
        wald.statistic = wald_statistic,
        wald.distribution = wald_distribution,
        wald.df1 = wald_df1,
        wald.df2 = wald_df2,
        wald.p.value = wald_p_value,
        pseudo.r.squared = .econ_multinomial_pseudo_r2(object),
        extra.significance = extra.significance
    )

    class(out) = "summary.econ_multinomial"

    out
}



#' @export
confint.econ_multinomial = function(
    object, parm = NULL, level = 0.95, ...
) {

    if (!is.numeric(level) ||
        length(level) != 1L ||
        !is.finite(level) ||
        level <= 0 ||
        level >= 1) {
        stop(
            "'level' must be strictly between 0 and 1.",
            call. = FALSE
        )
    }

    C = stats::coef(object)

    if (is.null(dim(C))) {
        C = matrix(C, nrow = 1L)
    }

    V = .econ_multinomial_beta_vcov(object)

    estimates = as.vector(t(C))
    parameter_names = colnames(V)

    if (is.null(parameter_names) ||
        length(estimates) != ncol(V) ||
        nrow(V) != ncol(V) ||
        !identical(rownames(V), parameter_names)) {
        stop(
            "Coefficient and covariance dimensions or names differ.",
            call. = FALSE
        )
    }

    names(estimates) = parameter_names

    if (is.null(parm)) {
        selected = seq_along(estimates)
    } else {

        if (!is.character(parm) ||
            anyNA(parm) ||
            anyDuplicated(parm) ||
            !all(parm %in% parameter_names)) {
            stop(
                "'parm' must contain valid, unique coefficient names.",
                call. = FALSE
            )
        }

        selected = match(parm, parameter_names)
    }

    se = sqrt(diag(V))

    clustered = object$vcov_type %in% c(
        "CR0", "CR1", "CR1p", "CR1S", "stata"
    )

    critical = if (clustered) {
        G = length(unique(object$cluster))
        stats::qt(1 - (1 - level) / 2, df = G - 1L)
    } else {
        stats::qnorm(1 - (1 - level) / 2)
    }

    result = cbind(
        estimates[selected] - critical * se[selected],
        estimates[selected] + critical * se[selected]
    )

    rownames(result) = parameter_names[selected]

    colnames(result) = paste0(
        format(
            100 * c((1 - level) / 2, 1 - (1 - level) / 2),
            trim = TRUE
        ),
        " %"
    )

    result
}


###### 4.2. Print Multinomial Summary ####

#' @export
print.summary.econ_multinomial = function(x, ...) {

    ###### 1. Print call and model information ####

    cat("\nCall:\n")
    print(x$call)

    cat(
        "\nMultinomial response model: ", x$link,
        "\nCovariance type: ", x$vcov_type,
        "\nObservations: ", x$nobs,
        "\n",
        sep = ""
    )

    if (is.finite(x$df)) {
        cat(
            "Inference: t distribution (df = ",
            x$df, ")\n",
            sep = ""
        )
    } else {
        cat("Inference: standard normal distribution\n")
    }

    ###### 2. Print coefficient inference ####

    for (category in rownames(x$coefficients)) {

        coef_table = cbind(
            Estimate = x$coefficients[category, ],
            "Std. Error" = x$standard.errors[category, ],
            statistic = x$statistic[category, ],
            p.value = x$p.value[category, ]
        )

        colnames(coef_table)[3L] =
            if (is.finite(x$df)) "t value" else "z value"

        colnames(coef_table)[4L] =
            if (is.finite(x$df)) "Pr(>|t|)" else "Pr(>|z|)"

        p_values = x$p.value[category, ]

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

        formatted = matrix(
            sprintf("%.5f", coef_table),
            nrow = nrow(coef_table),
            dimnames = dimnames(coef_table)
        )

        formatted = cbind(formatted, " " = stars)

        cat("\nOutcome:", category, "\n")

        print(
            formatted,
            quote = FALSE,
            right = TRUE
        )
    }

    ###### 3. Print significance legend ####

    if (isTRUE(x$extra.significance)) {

        cat(
            "\nSignificance: * p < 0.1, ** p < 0.05, ",
            "*** p < 0.01, + p < 0.001\n",
            sep = ""
        )

    } else {

        cat(
            "\nSignificance: * p < 0.1, ** p < 0.05, ",
            "*** p < 0.01\n",
            sep = ""
        )
    }

    ###### 4. Print joint Wald test ####

    cat("\nJoint Wald test (all slopes = 0):\n")

    if (is.finite(x$wald.statistic)) {

        if (x$wald.distribution == "F") {

            cat(
                "F(", x$wald.df1, ", ", x$wald.df2, ") = ",
                format(x$wald.statistic, digits = 5),
                ", p-value = ",
                format.pval(x$wald.p.value),
                "\n",
                sep = ""
            )

        } else {

            cat(
                "Chi-squared(", x$wald.df1, ") = ",
                format(x$wald.statistic, digits = 5),
                ", p-value = ",
                format.pval(x$wald.p.value),
                "\n",
                sep = ""
            )
        }

    } else {

        cat("Not available\n")
    }

    ###### 5. Print model fit ####

    cat(
        "\nLog-likelihood: ",
        sprintf("%.5f", as.numeric(x$logLik)),
        "\nMcFadden pseudo-R-squared: ",
        sprintf("%.5f", x$pseudo.r.squared),
        "\n",
        sep = ""
    )

    invisible(x)
}



###### 5. Coefficients, Predictions, and Model Information ####


#' @export
coef.econ_multinomial = function(object, ...) {

    if (identical(object$link, "probit")) {
        return(object$coefficients)
    }

    NextMethod()
}


#' @export
predict.econ_multinomial = function(
    object,
    newdata = NULL,
    type = c("class", "probs"),
    ...
) {

    type = match.arg(type)

    if (!identical(object$link, "probit")) {
        return(NextMethod())
    }

    if (is.null(newdata)) {

        P = object$fitted.values

    } else {

        terms_object = stats::delete.response(object$terms)

        mf_new = stats::model.frame(
            terms_object,
            data = newdata,
            na.action = stats::na.pass,
            xlev = object$xlevels
        )

        X_new = stats::model.matrix(
            terms_object,
            data = mf_new,
            contrasts.arg = object$contrasts
        )

        if (!identical(
            colnames(X_new),
            colnames(object$coefficients)
        )) {
            stop(
                "Prediction model matrix does not match fitted coefficients.",
                call. = FALSE
            )
        }

        complete = stats::complete.cases(X_new)

        P = matrix(
            NA_real_,
            nrow = nrow(X_new),
            ncol = length(object$lev),
            dimnames = list(NULL, object$lev)
        )

        if (any(complete)) {
            
            eta = cbind(
                0,
                X_new[complete, , drop = FALSE] %*%
                    t(object$coefficients)
            )
            
            colnames(eta) = object$lev
            
            P[complete, ] =
                .econ_multinomial_probit_probs_general(        
            eta = eta,
            Omega = object$Omega,
            normalize = TRUE
        )
        }
    }

    if (type == "probs") {
        return(P)
    }

    selected = rep(NA_integer_, nrow(P))
    complete = stats::complete.cases(P)

    if (any(complete)) {
        selected[complete] = max.col(
            P[complete, , drop = FALSE],
            ties.method = "first"
        )
    }

    factor(
        object$lev[selected],
        levels = object$lev
    )
}


#' @export
fitted.econ_multinomial = function(object, ...) {

    if (identical(object$link, "probit")) {
        return(object$fitted.values)
    }

    NextMethod()
}


#' @export
logLik.econ_multinomial = function(object, ...) {

    if (!identical(object$link, "probit")) {
        return(NextMethod())
    }

    structure(
        object$logLik,
        class = "logLik",
        df = length(object$coefficients) + length(object$covariance_parameters),
        nobs = length(object$used)
    )
}


#' @export
nobs.econ_multinomial = function(object, ...) length(object$used)

#' @export



###### 6. Average Marginal Effects ####

#' Average Marginal Effects for Multinomial Logit
#'
#' Computes average marginal effects for continuous predictors and
#' average discrete probability changes for factor predictors.
#' Standard errors use the delta method and the covariance estimator
#' selected when fitting the model.
#'
#' @param object An econ_multinomial model.
#' @param variables Optional predictor names. By default, all supported
#'   numeric and factor predictors are included.
#' @param level Confidence level strictly between 0 and 1.
#' @param weighted Logical. If TRUE, average using estimation weights.
#' @param ... Unused.
#' @return A data frame containing outcome, variable, contrast,
#'   AME, standard error, statistic, p-value and confidence intervals.
#' @export
marginal_effects.econ_multinomial = function(
    object,
    variables = NULL,
    level = 0.95,
    weighted = FALSE,
    ...
) {

    ###### 1. Validate inputs ####

    if (!inherits(object, "econ_multinomial")) {
        stop("'object' must be an econ_multinomial model.",
             call. = FALSE)
    }

    if (!is.numeric(level) || length(level) != 1L ||
        !is.finite(level) || level <= 0 || level >= 1) {
        stop("'level' must be strictly between 0 and 1.",
             call. = FALSE)
    }

    if (!is.logical(weighted) ||
        length(weighted) != 1L || is.na(weighted)) {
        stop("'weighted' must be TRUE or FALSE.",
             call. = FALSE)
    }

    mf = object$econ_model
    terms_object = stats::delete.response(object$terms)
    predictors = all.vars(terms_object)

    supported = vapply(
        predictors,
        function(v) is.numeric(mf[[v]]) || is.factor(mf[[v]]),
        logical(1)
    )

    if (is.null(variables)) {
        variables = predictors[supported]
    }

    if (!is.character(variables) ||
        length(variables) == 0L ||
        anyNA(variables) ||
        anyDuplicated(variables) ||
        !all(variables %in% predictors) ||
        !all(variables %in% predictors[supported])) {
        stop("'variables' must contain valid, unique numeric or factor predictors.",
             call. = FALSE)
    }

    term_labels = attr(terms_object, "term.labels")

    if (any(grepl(":", term_labels, fixed = TRUE))) {
        stop("Interactions are not yet supported for marginal effects.",
             call. = FALSE)
    }

    if (!all(variables %in% term_labels)) {
        stop("Marginal effects require additive, untransformed predictors.",
             call. = FALSE)
    }

    ###### 2. Extract model components ####

    X = stats::model.matrix(
        terms_object,
        data = mf,
        contrasts.arg = object$contrasts
    )

    C = stats::coef(object)

    if (is.null(dim(C))) {
        C = matrix(
            C,
            nrow = 1L,
            dimnames = list(
                object$lev[2L],
                names(C)
            )
        )
    }

    outcomes = c(object$lev[1L], rownames(C))

    B = matrix(
        0,
        nrow = length(outcomes),
        ncol = ncol(X),
        dimnames = list(outcomes, colnames(X))
    )

    B[-1L, ] = C[, colnames(X), drop = FALSE]

    V = stats::vcov(object)

    beta_names = as.vector(t(outer(
        rownames(C),
        colnames(C),
        paste,
        sep = ":"
    )))

    beta = as.vector(t(C))
    names(beta) = beta_names

    gamma = if (identical(object$link, "probit") &&
                identical(object$error_covariance, "estimated")) {
        object$covariance_parameters
    } else {
        numeric(0)
    }
    theta = c(beta, gamma)

    if (!setequal(names(theta), colnames(V)) ||
        !identical(rownames(V), colnames(V))) {
        stop("Parameter and covariance names differ.", call. = FALSE)
    }
    theta = theta[colnames(V)]

    ###### 3. Set averaging weights ####

    if (weighted) {

        averaging_weights = object$econ_weights

        if (is.null(averaging_weights) ||
            length(averaging_weights) != nrow(X) ||
            anyNA(averaging_weights) ||
            any(!is.finite(averaging_weights)) ||
            any(averaging_weights < 0) ||
            sum(averaging_weights) <= 0) {
            stop("Valid estimation weights are required.",
                 call. = FALSE)
        }

        averaging_weights = averaging_weights /
            sum(averaging_weights)

    } else {

        averaging_weights = rep(1 / nrow(X), nrow(X))
    }

    ###### 4. Construct factor contrasts ####

    contrast_data = list()

    for (v in variables) {

        if (!is.factor(mf[[v]])) {
            contrast_data[[v]] = list(
                list(
                    label = NA_character_,
                    X_high = NULL,
                    X_low = NULL
                )
            )

            next
        }

        lev = levels(mf[[v]])

        if (length(lev) < 2L) {
            stop("Factor predictors require at least two levels.",
                 call. = FALSE)
        }

        reference = lev[1L]

        contrast_data[[v]] = lapply(
            lev[-1L],
            function(target) {

                high = mf
                low = mf

                high[[v]] = factor(
                    rep(target, nrow(mf)),
                    levels = lev,
                    ordered = is.ordered(mf[[v]])
                )

                low[[v]] = factor(
                    rep(reference, nrow(mf)),
                    levels = lev,
                    ordered = is.ordered(mf[[v]])
                )

                list(
                    label = paste0(target, " - ", reference),
                    X_high = stats::model.matrix(
                        terms_object,
                        data = high,
                        contrasts.arg = object$contrasts
                    ),
                    X_low = stats::model.matrix(
                        terms_object,
                        data = low,
                        contrasts.arg = object$contrasts
                    )
                )
            }
        )
    }

    ###### 5. Build result labels ####

    labels = do.call(
        rbind,
        lapply(variables, function(v) {

            data.frame(
                outcome = rep(
                    outcomes,
                    times = length(contrast_data[[v]])
                ),
                variable = v,
                contrast = rep(
                    vapply(
                        contrast_data[[v]],
                        function(z) z$label,
                        character(1)
                    ),
                    each = length(outcomes)
                ),
                stringsAsFactors = FALSE
            )
        })
    )

    rownames(labels) = NULL

###### 6. Compute effects from coefficients ####

probabilities = function(
    X_new,
    B_current,
    theta_current
) {

    eta = X_new %*% t(B_current)

    if (identical(object$link, "probit")) {

        return(
            .econ_multinomial_probit_probs_general(
                eta = eta,
                normalize = TRUE,
                Omega = if (length(gamma) > 0L) {
                    .econ_multinomial_probit_omega(
                        theta_current[names(gamma)],
                        length(outcomes) - 1L
                    )
                } else {
                    object$Omega
                }
            )
        )
    }

    eta = sweep(
        eta,
        1L,
        apply(eta, 1L, max),
        "-"
    )

    exp_eta = exp(eta)

    exp_eta / rowSums(exp_eta)
}


    calculate_ame = function(theta_current) {

        b = theta_current[beta_names]

        C_current = matrix(
            b[beta_names],
            nrow = nrow(C),
            byrow = TRUE,
            dimnames = dimnames(C)
        )

        B_current = B
        B_current[-1L, ] = C_current[
            ,
            colnames(X),
            drop = FALSE
        ]

        P = probabilities(X, B_current, theta_current)

        results = list()

        for (v in variables) {

            if (!is.factor(mf[[v]])) {

                if (!(v %in% colnames(X))) {
                    stop("Unsupported numeric predictor transformation.",
                         call. = FALSE)
                }
              
              
                if (identical(object$link, "probit")) {

                    h = 1e-5

                    X_plus = X
                    X_minus = X

                    X_plus[, v] = X_plus[, v] + h
                    X_minus[, v] = X_minus[, v] - h

                    P_plus = probabilities(X_plus, B_current, theta_current)
                    P_minus = probabilities(X_minus, B_current, theta_current)

                    derivatives = (
                        P_plus - P_minus
                    ) / (2 * h)

                    results[[length(results) + 1L]] =
                        colSums(derivatives * averaging_weights)

                    next
                }

                slopes = B_current[, v]

                mean_slope = as.vector(P %*% slopes)

                derivatives = P * sweep(
                    matrix(
                        slopes,
                        nrow = nrow(P),
                        ncol = ncol(P),
                        byrow = TRUE
                    ),
                    1L,
                    mean_slope,
                    "-"
                )

                results[[length(results) + 1L]] =
                    colSums(derivatives * averaging_weights)

            } else {

                for (comparison in contrast_data[[v]]) {

                    P_high = probabilities(
                        comparison$X_high,
                        B_current,
                        theta_current
                    )

                    P_low = probabilities(
                        comparison$X_low,
                        B_current,
                        theta_current
                    )

                    difference = P_high - P_low

                    results[[length(results) + 1L]] =
                        colSums(difference * averaging_weights)
                }
            }
        }

        unlist(results, use.names = FALSE)
    }

    ###### 7. Calculate average marginal effects ####

    estimates = calculate_ame(theta)

    ###### 8. Numerical delta-method Jacobian ####

    gradient = matrix(
        0,
        nrow = length(estimates),
        ncol = length(theta)
    )

    for (j in seq_along(theta)) {

        h = 1e-5 * max(1, abs(theta[j]))

        theta_plus = theta
        theta_minus = theta

        theta_plus[j] = theta_plus[j] + h
        theta_minus[j] = theta_minus[j] - h

        gradient[, j] = (
            calculate_ame(theta_plus) -
                calculate_ame(theta_minus)
        ) / (2 * h)
    }

    ###### 9. Delta-method inference ####

    V_ame = gradient %*% V %*% t(gradient)

    se = sqrt(pmax(0, diag(V_ame)))

    clustered = object$vcov_type %in% c(
        "CR0", "CR1", "CR1p", "CR1S", "stata"
    )

    df = if (clustered) {
        length(unique(object$cluster)) - 1L
    } else {
        Inf
    }

    statistic = estimates / se

    p_value = if (is.finite(df)) {
        2 * stats::pt(-abs(statistic), df = df)
    } else {
        2 * stats::pnorm(-abs(statistic))
    }

    critical = if (is.finite(df)) {
        stats::qt(1 - (1 - level) / 2, df = df)
    } else {
        stats::qnorm(1 - (1 - level) / 2)
    }

    ###### 10. Return results ####

    data.frame(
        outcome = labels$outcome,
        variable = labels$variable,
        contrast = labels$contrast,
        AME = estimates,
        std.error = se,
        statistic = statistic,
        p.value = p_value,
        conf.low = estimates - critical * se,
        conf.high = estimates + critical * se,
        row.names = NULL
    )
}





###### 7. Multinomial Augment Method ####

#' Augment Multinomial Logit Models
#'
#' Returns predicted outcome categories and category probabilities.
#' Optionally computes delta-method standard errors using the
#' selected model covariance estimator.
#'
#' @param x An econ_multinomial model.
#' @param data Optional estimation data.
#' @param newdata Optional data for out-of-sample predictions.
#' @param se_fit Logical. Whether to calculate prediction standard errors.
#' @param ... Unused.
#'
#' @return A tibble containing data, predictions, and optionally
#'   standard errors for predicted category probabilities.
#' @export
augment.econ_multinomial = function(
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

    if (!is.null(newdata)) {

        prediction_data = as.data.frame(newdata)
        result = prediction_data

    } else {

        prediction_data = as.data.frame(x$econ_model)

        result = if (is.null(data)) {
            prediction_data
        } else {
            as.data.frame(data)
        }

        if (nrow(result) != nrow(prediction_data)) {
            stop(
                "'data' must contain the estimation observations.",
                call. = FALSE
            )
        }

        internal = c(
            "(weights)",
            ".econ_multi_row_id",
            ".econ_multi_weights",
            ".econ_multi_subset"
        )

        result = result[
            ,
            !names(result) %in% internal,
            drop = FALSE
        ]
    }

    ###### 3. Predict categories and probabilities ####

    probabilities = as.matrix(
        stats::predict(
            x,
            newdata = prediction_data,
            type = "probs"
        )
    )

    categories = colnames(probabilities)

    if (is.null(categories)) {
        stop(
            "Predicted probabilities must have category names.",
            call. = FALSE
        )
    }

    result$.fitted = stats::predict(
        x,
        newdata = prediction_data,
        type = "class"
    )

    probability_names = paste0(
        ".fitted_",
        make.names(categories, unique = TRUE)
    )

    for (j in seq_along(categories)) {
        result[[probability_names[j]]] =
            unname(probabilities[, j])
    }


    ###### 4. Delta-method prediction standard errors ####

    if (isTRUE(se_fit) && identical(x$link, "probit")) {

        se_prob = .econ_multinomial_probit_prediction_se(
            object = x,
            newdata = prediction_data
        )

        se_names = paste0(
            ".se.fit_",
            make.names(categories, unique = TRUE)
        )

        for (j in seq_along(categories)) {
            result[[se_names[j]]] = unname(se_prob[, j])
        }

        se_fit = FALSE
    }

    if (isTRUE(se_fit)) {


        ###### 4.1. Construct model matrix ####

        terms_object = stats::delete.response(
            stats::terms(x)
        )

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

        ###### 4.2. Extract coefficients and covariance ####

        C = stats::coef(x)

        if (is.null(dim(C))) {
            C = matrix(
                C,
                nrow = 1L,
                dimnames = list(
                    categories[2L],
                    names(C)
                )
            )
        }

        X = X[
            ,
            colnames(C),
            drop = FALSE
        ]

        nonbase = rownames(C)

        parameter_names = unlist(
            lapply(nonbase, function(category) {
                paste0(category, ":", colnames(X))
            }),
            use.names = FALSE
        )

        V = stats::vcov(x)

        if (!setequal(parameter_names, colnames(V))) {
            stop(
                "Coefficient and covariance names differ.",
                call. = FALSE
            )
        }

        V = V[
            parameter_names,
            parameter_names,
            drop = FALSE
        ]

        ###### 4.3. Calculate category probability gradients ####

        se_prob = matrix(
            NA_real_,
            nrow = nrow(X),
            ncol = length(categories),
            dimnames = list(NULL, categories)
        )

        for (j in seq_along(categories)) {

            category = categories[j]
            p_j = probabilities[, j]

            gradient = do.call(
                cbind,
                lapply(nonbase, function(r) {

                    multiplier = p_j * (
                        as.numeric(category == r) -
                            probabilities[, r]
                    )

                    X * multiplier
                })
            )

            ###### 4.4. Apply delta method ####

            variance = rowSums(
                (gradient %*% V) * gradient
            )

            se_prob[, j] = sqrt(
                pmax(0, variance)
            )
        }

        ###### 4.5. Add standard-error columns ####

        se_names = paste0(
            ".se.fit_",
            make.names(categories, unique = TRUE)
        )

        for (j in seq_along(categories)) {
            result[[se_names[j]]] =
                unname(se_prob[, j])
        }
    }

    ###### 5. Return augmented data ####

    tibble::as_tibble(result)
}




###### 8. McFadden Pseudo-R-Squared ####

.econ_multinomial_pseudo_r2 = function(object) {

    ###### 1. Extract outcome and weights ####

    y = stats::model.response(object$econ_model)
    w = object$econ_weights

    ###### 2. Calculate weighted category counts ####

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

    ###### 3. Calculate null log-likelihood ####

    probabilities = counts / sum(counts)

    null_loglik = sum(
        counts * log(probabilities)
    )

    ###### 4. Calculate McFadden pseudo-R-squared ####

    fitted_loglik = as.numeric(stats::logLik(object))

    if (!is.finite(fitted_loglik) ||
        !is.finite(null_loglik) ||
        null_loglik >= 0) {
        return(NA_real_)
    }

    1 - fitted_loglik / null_loglik
}


###### 8.1. Multinomial Model Glance Method ####

#' @export
glance.econ_multinomial = function(x, ...) {

  
    if (identical(x$link, "probit")) {

        result = tibble::tibble(
            edf = length(x$coefficients) + length(x$covariance_parameters),
            logLik = as.numeric(stats::logLik(x)),
            AIC = stats::AIC(x),
            BIC = stats::BIC(x),
            deviance = -2 * as.numeric(stats::logLik(x)),
            nobs = stats::nobs(x),
            pseudo.r.squared = .econ_multinomial_pseudo_r2(x)
        )

        return(result)
    }


    ###### 1. Obtain standard multinomial model statistics ####

    result = utils::getS3method(
        "glance",
        "multinom"
    )(x, ...)

    ###### 2. Add McFadden pseudo-R-squared ####

    result$pseudo.r.squared =
        .econ_multinomial_pseudo_r2(x)

    ###### 3. Return model statistics ####

    result
}



###### 9. Multinomial Tidy Method ####

#' Tidy a Multinomial Logit Model
#'
#' Returns coefficient estimates and inference statistics using
#' the covariance estimator selected when fitting the model.
#'
#' @param x An econ_multinomial model.
#' @param conf.int Whether to include confidence intervals.
#' @param conf.level Confidence level.
#' @param ... Unused.
#'
#' @return A tibble containing coefficient estimates and inference statistics.
#' @export
tidy.econ_multinomial = function(
    x,
    conf.int = FALSE,
    conf.level = 0.95,
    ...
) {

    ###### 1. Extract EconR model summary ####

    s = summary(x)

    ###### 2. Construct coefficient table ####

    categories = rownames(s$coefficients)
    terms = colnames(s$coefficients)

    result = data.frame(
        y.level = rep(categories, each = length(terms)),
        term = rep(terms, times = length(categories)),
        estimate = as.vector(t(s$coefficients)),
        std.error = as.vector(t(s$standard.errors)),
        statistic = as.vector(t(s$statistic)),
        p.value = as.vector(t(s$p.value)),
        row.names = NULL
    )

    ###### 3. Add confidence intervals ####

    if (isTRUE(conf.int)) {

        if (!is.numeric(conf.level) ||
            length(conf.level) != 1L ||
            is.na(conf.level) ||
            !is.finite(conf.level) ||
            conf.level <= 0 ||
            conf.level >= 1) {

            stop(
                "'conf.level' must be strictly between 0 and 1.",
                call. = FALSE
            )
        }

        critical = if (is.finite(s$df)) {
            stats::qt(
                1 - (1 - conf.level) / 2,
                df = s$df
            )
        } else {
            stats::qnorm(
                1 - (1 - conf.level) / 2
            )
        }

        result$conf.low =
            result$estimate - critical * result$std.error

        result$conf.high =
            result$estimate + critical * result$std.error
    }

    ###### 4. Return tidy results ####

    tibble::as_tibble(result)
}



###### 10. Multinomial Probit Probability Engine ####

.econ_multinomial_probit_probs = function(eta, Omega) {

    ###### 1. Validate inputs ####

    eta = as.matrix(eta)

    if (!is.numeric(eta) ||
        ncol(eta) != 3L ||
        any(!is.finite(eta))) {
        stop(
            "'eta' must be a finite numeric matrix with three columns.",
            call. = FALSE
        )
    }

    if (!is.matrix(Omega) ||
        !is.numeric(Omega) ||
        !identical(dim(Omega), c(2L, 2L)) ||
        any(!is.finite(Omega)) ||
        !isTRUE(all.equal(Omega, t(Omega))) ||
        any(eigen(Omega, symmetric = TRUE)$values <= 0)) {
        stop(
            "'Omega' must be a symmetric positive-definite 2 x 2 matrix.",
            call. = FALSE
        )
    }

    ###### 2. Construct utility-error covariance ####

    Sigma = matrix(0, 3L, 3L)
    Sigma[-1L, -1L] = Omega

    probabilities = matrix(
        NA_real_,
        nrow = nrow(eta),
        ncol = 3L,
        dimnames = list(NULL, colnames(eta))
    )

    ###### 3. Evaluate choice probabilities ####

    for (j in seq_len(3L)) {

        others = setdiff(seq_len(3L), j)

        D = matrix(0, 2L, 3L)

        for (k in seq_along(others)) {
            D[k, others[k]] = 1
            D[k, j] = -1
        }

        Sigma_j = D %*% Sigma %*% t(D)

        for (i in seq_len(nrow(eta))) {

            upper = eta[i, j] - eta[i, others]

            probabilities[i, j] = as.numeric(
                mvtnorm::pmvnorm(
                    lower = rep(-Inf, 2L),
                    upper = upper,
                    mean = rep(0, 2L),
                    sigma = Sigma_j,
                    algorithm = mvtnorm::Miwa()
                )
            )
        }
    }

    probabilities
}



###### 10.1. Restricted Multinomial Probit Maximum Likelihood ####

.econ_multinomial_probit_fit = function(
    X,
    y,
    weights = NULL,
    maxit = 500
) {

    ###### 1. Validate inputs ####

    X = as.matrix(X)

    if (!is.numeric(X) ||
        nrow(X) == 0L ||
        ncol(X) == 0L ||
        any(!is.finite(X))) {
        stop("'X' must be a finite numeric model matrix.", call. = FALSE)
    }

    if (!is.factor(y) ||
        is.ordered(y) ||
        nlevels(y) != 3L ||
        length(y) != nrow(X) ||
        anyNA(y) ||
        any(tabulate(as.integer(y), nbins = 3L) == 0L)) {
        stop(
            "'y' must be an unordered three-category factor with all categories observed.",
            call. = FALSE
        )
    }

    if (is.null(weights)) {
        weights = rep(1, nrow(X))
    }

    if (!is.numeric(weights) ||
        length(weights) != nrow(X) ||
        anyNA(weights) ||
        any(!is.finite(weights)) ||
        any(weights < 0) ||
        sum(weights) <= 0) {
        stop("Invalid observation weights.", call. = FALSE)
    }

    ###### 2. Set identification restrictions ####

    categories = levels(y)
    nonbase = categories[-1L]

    Omega = matrix(
        c(1, 0.5, 0.5, 1),
        nrow = 2L
    )

    K = ncol(X)
    J = 3L

    ###### 3. Define parameter mapping ####

    make_coefficients = function(beta) {

        matrix(
            beta,
            nrow = J - 1L,
            ncol = K,
            byrow = TRUE,
            dimnames = list(nonbase, colnames(X))
        )
    }

    probability_function = function(beta) {

        B = make_coefficients(beta)

        eta = cbind(
            0,
            X %*% t(B)
        )

        colnames(eta) = categories

        .econ_multinomial_probit_probs(
            eta = eta,
            Omega = Omega
        )
    }

    ###### 4. Define weighted negative log-likelihood ####

    outcome_index = cbind(
        seq_len(nrow(X)),
        as.integer(y)
    )

    negative_loglik = function(beta) {

        P = probability_function(beta)

        observed_probabilities = P[outcome_index]

        if (any(!is.finite(observed_probabilities)) ||
            any(observed_probabilities < 0)) {
            return(.Machine$double.xmax / 100)
        }

        observed_probabilities = pmax(
            observed_probabilities,
            .Machine$double.xmin
        )

        -sum(weights * log(observed_probabilities))
    }

    ###### 5. Estimate model ####

    initial = rep(0, (J - 1L) * K)

    optimization = stats::optim(
        par = initial,
        fn = negative_loglik,
        method = "BFGS",
        control = list(
            maxit = maxit,
            reltol = 1e-9
        )
    )

    if (optimization$convergence != 0L) {
        warning(
            "Multinomial probit optimization did not converge.",
            call. = FALSE
        )
    }

        ###### 6. Calculate Hessian and classical covariance ####

    hessian = stats::optimHess(
        par = optimization$par,
        fn = negative_loglik
    )

    hessian = (hessian + t(hessian)) / 2

    parameter_names = unlist(
        lapply(nonbase, function(category) {
            paste0(category, ":", colnames(X))
        }),
        use.names = FALSE
    )

    dimnames(hessian) = list(
        parameter_names,
        parameter_names
    )

    if (any(!is.finite(hessian)) ||
        min(eigen(hessian, symmetric = TRUE)$values) <= 0) {
        stop(
            "Multinomial probit Hessian is not positive definite.",
            call. = FALSE
        )
    }

    V = solve(hessian)

    dimnames(V) = dimnames(hessian)

    
    ###### 7. Construct estimation results ####

    coefficients = make_coefficients(
        optimization$par
    )

    fitted_probabilities = probability_function(
        optimization$par
    )

    list(
        coefficients = coefficients,
        vcov = V,
        hessian = hessian,
        logLik = -optimization$value,
        fitted.values = fitted_probabilities,
        weights = weights,
        Omega = Omega,
        convergence = optimization$convergence,
        counts = optimization$counts
    )
}




###### 11. Generalized Multinomial Probit ####


###### 11.1. Generalized Choice Probabilities ####

.econ_multinomial_probit_probs_general = function(
    eta,
    Omega,
    observed = NULL,
    normalize = FALSE
) {

    eta = as.matrix(eta)
    Omega = as.matrix(Omega)

    J = ncol(eta)
    M = J - 1L

    if (!is.numeric(eta) ||
        J < 3L ||
        nrow(eta) < 1L ||
        any(!is.finite(eta))) {
        stop("Invalid systematic utility matrix.", call. = FALSE)
    }

    if (J > 20L) {
        stop(
            "This deterministic Miwa implementation supports at most 20 alternatives.",
            call. = FALSE
        )
    }

    if (!is.numeric(Omega) ||
        !identical(dim(Omega), c(M, M)) ||
        any(!is.finite(Omega)) ||
        !isTRUE(all.equal(Omega, t(Omega))) ||
        min(eigen(Omega, symmetric = TRUE)$values) <= 0) {
        stop(
            "'Omega' must be a positive-definite (J - 1) x (J - 1) matrix.",
            call. = FALSE
        )
    }

    if (!is.null(observed)) {

        if (length(observed) != nrow(eta) ||
            anyNA(observed) ||
            any(!(observed %in% seq_len(J)))) {
            stop(
                "Invalid observed choice indices.",
                call. = FALSE
            )
        }

        observed = as.integer(observed)
    }

    if (!is.logical(normalize) ||
        length(normalize) != 1L ||
        is.na(normalize) ||
        (normalize && !is.null(observed))) {
        stop(
            "Invalid probability normalization request.",
            call. = FALSE
        )
    }

    ###### 1. Construct normalized utility covariance ####

    Sigma = matrix(0, J, J)
    Sigma[-1L, -1L] = Omega

    P = matrix(
        NA_real_,
        nrow = nrow(eta),
        ncol = J,
        dimnames = list(NULL, colnames(eta))
    )

    if (!is.null(observed)) {
        result = numeric(nrow(eta))
    }

    ###### 2. Calculate category probabilities ####

    for (j in seq_len(J)) {

        others = setdiff(seq_len(J), j)

        D = matrix(0, M, J)

        for (k in seq_along(others)) {

            D[k, others[k]] = 1
            D[k, j] = -1
        }

        Sigma_j = D %*% Sigma %*% t(D)

        rows = if (is.null(observed)) {
            seq_len(nrow(eta))
        } else {
            which(observed == j)
        }

        for (i in rows) {

            upper = eta[i, j] - eta[i, others]

            value = as.numeric(
                mvtnorm::pmvnorm(
                    lower = rep(-Inf, M),
                    upper = upper,
                    mean = rep(0, M),
                    sigma = Sigma_j,
                    algorithm = mvtnorm::Miwa()
                )
            )

            if (is.null(observed)) {
                P[i, j] = value
            } else {
                result[i] = value
            }
        }
    }

    ###### 3. Return observed or complete probabilities ####

    if (!is.null(observed)) {
        return(result)
    }

    if (normalize) {

        totals = rowSums(P)

        if (any(!is.finite(totals)) ||
            any(totals <= 0)) {
            stop(
                "Invalid multinomial probit probability totals.",
                call. = FALSE
            )
        }

        P = P / totals
    }

    P
}


###### 11.2. Identified Multinomial Probit Error Covariance ####

.econ_multinomial_probit_omega = function(
    gamma,
    M
) {

    ###### 1. Validate covariance dimensions ####

    if (!is.numeric(M) ||
        length(M) != 1L ||
        !is.finite(M) ||
        M != as.integer(M) ||
        M < 2L) {
        stop(
            "'M' must be an integer greater than or equal to two.",
            call. = FALSE
        )
    }

    M = as.integer(M)

    q = M * (M + 1L) / 2L - 1L

    if (!is.numeric(gamma) ||
        length(gamma) != q ||
        any(!is.finite(gamma))) {
        stop(
            "Invalid multinomial probit covariance parameters.",
            call. = FALSE
        )
    }

    ###### 2. Construct normalized Cholesky factor ####

    L = matrix(
        0,
        nrow = M,
        ncol = M
    )

    L[1L, 1L] = 1

    position = 1L

    for (i in seq.int(2L, M)) {

        for (j in seq_len(i)) {

            if (i == j) {

                L[i, j] = exp(gamma[position])

            } else {

                L[i, j] = gamma[position]
            }

            position = position + 1L
        }
    }

    ###### 3. Construct covariance ####

    Omega = tcrossprod(L)

    if (any(!is.finite(Omega)) ||
        any(!is.finite(L)) ||
        any(diag(L) <= 0)) {
        stop(
            "Covariance parameters produce an invalid matrix.",
            call. = FALSE
        )
    }

    Omega
}


###### 11.3. Multinomial Probit Covariance Parameter Names ####

.econ_multinomial_probit_cov_names = function(M) {

    names = character(0)

    for (i in seq.int(2L, M)) {

        for (j in seq_len(i)) {

            names = c(
                names,
                if (i == j) {
                    paste0("cov:logL[", i, ",", j, "]")
                } else {
                    paste0("cov:L[", i, ",", j, "]")
                }
            )
        }
    }

    names
}


###### 11.4. Generalized Multinomial Probit MLE ####

.econ_multinomial_probit_fit_general = function(
    X,
    y,
    weights = NULL,
    maxit = 500,
    error_covariance = c("fixed", "estimated")
) {

    error_covariance = match.arg(error_covariance)

    ###### 1. Validate inputs ####

    X = as.matrix(X)

    if (!is.numeric(X) ||
        nrow(X) < 1L ||
        ncol(X) < 1L ||
        any(!is.finite(X))) {
        stop(
            "'X' must be a finite numeric model matrix.",
            call. = FALSE
        )
    }

    if (!is.factor(y) ||
        is.ordered(y) ||
        nlevels(y) < 3L ||
        length(y) != nrow(X) ||
        anyNA(y) ||
        any(tabulate(
            as.integer(y),
            nbins = nlevels(y)
        ) == 0L)) {
        stop(
            "'y' must be an unordered factor with at least three observed categories.",
            call. = FALSE
        )
    }

    if (is.null(weights)) {
        weights = rep(1, nrow(X))
    }

    if (!is.numeric(weights) ||
        length(weights) != nrow(X) ||
        anyNA(weights) ||
        any(!is.finite(weights)) ||
        any(weights < 0) ||
        sum(weights) <= 0) {
        stop(
            "Invalid observation weights.",
            call. = FALSE
        )
    }

    ###### 2. Identification and parameter setup ####

    categories = levels(y)

    J = length(categories)
    M = J - 1L
    K = ncol(X)

    if (J > 20L) {
        stop(
            "This implementation supports at most 20 alternatives.",
            call. = FALSE
        )
    }

    if (is.null(colnames(X))) {
        colnames(X) = paste0("x", seq_len(K))
    }

    ###### 2.1. Fixed covariance reference specification ####

    # Independent normal utility errors with variance 0.5.
    # Utility differences relative to category 1 have
    # variance 1 and pairwise covariance 0.5.

    Omega_fixed = matrix(0.5, M, M)
    diag(Omega_fixed) = 1

    ###### 2.2. Utility coefficient parameter names ####

    beta_names = unlist(
        lapply(categories[-1L], function(category) {
            paste0(category, ":", colnames(X))
        }),
        use.names = FALSE
    )

    n_beta = M * K

    ###### 2.3. Covariance parameter setup ####

    if (error_covariance == "estimated") {

        cov_names = .econ_multinomial_probit_cov_names(M)

        L_initial = t(chol(Omega_fixed))

        gamma_initial = numeric(0)

        for (i in seq.int(2L, M)) {

            for (j in seq_len(i)) {

                gamma_initial = c(
                    gamma_initial,
                    if (i == j) {
                        log(L_initial[i, j])
                    } else {
                        L_initial[i, j]
                    }
                )
            }
        }

    } else {

        cov_names = character(0)
        gamma_initial = numeric(0)
    }

    parameter_names = c(
        beta_names,
        cov_names
    )

    ###### 2.4. Construct covariance for parameters ####

    get_omega = function(theta) {

        if (error_covariance == "fixed") {
            return(Omega_fixed)
        }

        gamma = theta[
            seq.int(n_beta + 1L, length(theta))
        ]

        .econ_multinomial_probit_omega(
            gamma = gamma,
            M = M
        )
    }

    ###### 3. Construct coefficients and probabilities ####

    make_coefficients = function(beta) {

        matrix(
            beta,
            nrow = M,
            ncol = K,
            byrow = TRUE,
            dimnames = list(
                categories[-1L],
                colnames(X)
            )
        )
    }

    probability_function = function(
        theta,
        observed = NULL,
        normalize = FALSE
    ) {

        beta = theta[seq_len(n_beta)]

        B = make_coefficients(beta)

        eta = cbind(
            0,
            X %*% t(B)
        )

        colnames(eta) = categories

        Omega = get_omega(theta)

        .econ_multinomial_probit_probs_general(
            eta = eta,
            Omega = Omega,
            observed = observed,
            normalize = normalize
        )
    }

    ###### 4. Weighted negative log-likelihood ####

    outcome_index = cbind(
        seq_len(nrow(X)),
        as.integer(y)
    )

    negative_loglik = function(theta) {

        observed = tryCatch(
            probability_function(
                theta,
                observed = as.integer(y)
            ),
            error = function(e) NULL
        )

        if (is.null(observed)) {
            return(.Machine$double.xmax / 100)
        }

        if (any(!is.finite(observed)) ||
            any(observed < 0)) {
            return(.Machine$double.xmax / 100)
        }

        observed = pmax(
            observed,
            .Machine$double.xmin
        )

        -sum(weights * log(observed))
    }

    ###### 5. Maximum-likelihood estimation ####

    initial = c(
        rep(0, n_beta),
        gamma_initial
    )

    names(initial) = parameter_names

    optimization = stats::optim(
        par = initial,
        fn = negative_loglik,
        method = "BFGS",
        control = list(
            maxit = maxit,
            reltol = 1e-9
        )
    )

    if (optimization$convergence != 0L) {
        warning(
            "Generalized multinomial probit optimization did not converge.",
            call. = FALSE
        )
    }

    ###### 6. Classical covariance ####

    hessian = stats::optimHess(
        par = optimization$par,
        fn = negative_loglik
    )

    hessian = (hessian + t(hessian)) / 2

    dimnames(hessian) = list(
        parameter_names,
        parameter_names
    )

    if (any(!is.finite(hessian)) ||
        min(eigen(
            hessian,
            symmetric = TRUE
        )$values) <= 0) {
        stop(
            "Multinomial probit Hessian is not positive definite.",
            call. = FALSE
        )
    }

    V = solve(hessian)

    dimnames(V) = dimnames(hessian)

    ###### 7. Return estimation results ####

    theta_hat = optimization$par

    beta_hat = theta_hat[seq_len(n_beta)]

    gamma_hat = if (length(cov_names) > 0L) {
        theta_hat[
            seq.int(n_beta + 1L, length(theta_hat))
        ]
    } else {
        numeric(0)
    }

    coefficients = make_coefficients(
        beta_hat
    )

    Omega = get_omega(theta_hat)

    fitted_probabilities = probability_function(
        theta_hat,
        normalize = TRUE
    )

    names(gamma_hat) = cov_names

    list(
        coefficients = coefficients,
        vcov = V,
        hessian = hessian,
        logLik = -optimization$value,
        fitted.values = fitted_probabilities,
        weights = weights,
        Omega = Omega,
        error_covariance = error_covariance,
        covariance_parameters = gamma_hat,
        parameter_names = parameter_names,
        convergence = optimization$convergence,
        counts = optimization$counts
    )
}


###### 12. Multinomial Probit Robust Covariance ####

###### 12.1. Observation-Level Score Contributions ####

.econ_multinomial_probit_scores = function(
    fit,
    X,
    y,
    weights = NULL,
    step = 1e-5
) {

    X = as.matrix(X)

    if (!is.numeric(X) ||
        any(!is.finite(X)) ||
        !is.factor(y) ||
        is.ordered(y) ||
        nrow(X) != length(y) ||
        anyNA(y) ||
        !identical(
            levels(y),
            colnames(fit$fitted.values)
        )) {
        stop(
            "Invalid or misaligned score inputs.",
            call. = FALSE
        )
    }

    n = nrow(X)
    K = ncol(X)
    M = nlevels(y) - 1L

    if (is.null(weights)) {
        weights = fit$weights
    }

    if (!is.numeric(weights) ||
        length(weights) != n ||
        anyNA(weights) ||
        any(!is.finite(weights)) ||
        any(weights < 0)) {
        stop(
            "Invalid score weights.",
            call. = FALSE
        )
    }

    if (!is.numeric(step) ||
        length(step) != 1L ||
        !is.finite(step) ||
        step <= 0) {
        stop(
            "'step' must be positive and finite.",
            call. = FALSE
        )
    }

    ###### 1. Joint parameter vector ####

    n_beta = M * K

    beta = as.vector(
        t(fit$coefficients)
    )

    estimated = identical(
        fit$error_covariance,
        "estimated"
    )

    gamma = if (estimated) {
        fit$covariance_parameters
    } else {
        numeric(0)
    }

    theta = c(beta, gamma)

    if (length(beta) != n_beta ||
        length(theta) != ncol(fit$vcov)) {
        stop(
            "Probit parameter and covariance dimensions do not match.",
            call. = FALSE
        )
    }

    names(theta) = colnames(fit$vcov)

    ###### 2. Individual log-likelihood ####

    individual_loglik = function(par) {

        B = matrix(
            par[seq_len(n_beta)],
            nrow = M,
            ncol = K,
            byrow = TRUE
        )

        eta = cbind(
            0,
            X %*% t(B)
        )

        colnames(eta) = levels(y)

        Omega = if (estimated) {
            .econ_multinomial_probit_omega(
                par[-seq_len(n_beta)],
                M
            )
        } else {
            fit$Omega
        }

        observed_prob = .econ_multinomial_probit_probs_general(
            eta = eta,
            Omega = Omega,
            observed = as.integer(y)
        )

        log(
            pmax(
                observed_prob,
                .Machine$double.xmin
            )
        )
    }

    ###### 3. Numerical score contributions ####

    S = matrix(
        NA_real_,
        nrow = n,
        ncol = length(theta),
        dimnames = list(
            NULL,
            names(theta)
        )
    )

    for (j in seq_along(theta)) {

        h = step * max(
            1,
            abs(theta[j])
        )

        plus = theta
        minus = theta

        plus[j] = plus[j] + h
        minus[j] = minus[j] - h

        S[, j] = weights *
            (
                individual_loglik(plus) -
                individual_loglik(minus)
            ) / (2 * h)
    }

    if (any(!is.finite(S))) {
        stop(
            "Nonfinite probit score contributions.",
            call. = FALSE
        )
    }

    S
}

###### 12.2. HC and Cluster-Robust Covariance ####

.econ_multinomial_probit_vcov = function(
    fit,
    X,
    y,
    vcov = c(
        "classical", "HC0", "HC1",
        "CR0", "CR1", "CR1p", "CR1S", "stata"
    ),
    cluster = NULL,
    weights = NULL,
    scores = NULL
) {

    vcov = match.arg(vcov)

    V = fit$vcov

    if (vcov == "classical") {
        return(V)
    }

    ###### 1. Obtain or reuse score contributions ####

    if (is.null(scores)) {

        S = .econ_multinomial_probit_scores(
            fit = fit,
            X = X,
            y = y,
            weights = weights
        )

    } else {

        S = as.matrix(scores)

        if (!is.null(weights)) {
            stop(
                "Supply either precomputed scores or weights, not both.",
                call. = FALSE
            )
        }

        if (!is.numeric(S) ||
            !identical(dim(S), c(nrow(X), ncol(V))) ||
            any(!is.finite(S))) {
            stop(
                "Invalid precomputed probit scores.",
                call. = FALSE
            )
        }
    }

    if (!identical(colnames(S), colnames(V))) {
        stop(
            "Score and covariance parameter names differ.",
            call. = FALSE
        )
    }

    n = nrow(S)
    k = ncol(S)

    ###### 2. Heteroskedasticity-robust covariance ####

    if (vcov %in% c("HC0", "HC1")) {

        M = crossprod(S)

        if (vcov == "HC1") {

            if (n <= k) {
                stop("HC1 requires n > k.", call. = FALSE)
            }

            M = M * n / (n - k)
        }

    ###### 3. Cluster-robust covariance ####

    } else {

        if (is.null(cluster) ||
            length(cluster) != n ||
            anyNA(cluster)) {
            stop(
                "Invalid or misaligned cluster variable.",
                call. = FALSE
            )
        }

        G = length(unique(cluster))

        if (G < 2L) {
            stop(
                "At least two clusters are required.",
                call. = FALSE
            )
        }

        U = rowsum(
            S,
            group = cluster,
            reorder = TRUE
        )

        M = crossprod(U)

        ###### 4. Existing EconR correction factors ####

        if (vcov == "CR1") {

            if (n <= k) {
                stop("CR1 requires n > k.", call. = FALSE)
            }

            M = M * (n - 1) / (n - k)

        } else if (vcov == "CR1p") {

            if (G <= k) {
                stop("CR1p requires G > k.", call. = FALSE)
            }

            M = M * G / (G - k)

        } else if (vcov %in% c("CR1S", "stata")) {

            if (n <= k) {
                stop(
                    "CR1S/stata requires n > k.",
                    call. = FALSE
                )
            }

            M = M * G * (n - 1) /
                ((G - 1) * (n - k))
        }
    }

    ###### 5. Sandwich covariance ####

    result = V %*% M %*% V

    dimnames(result) = dimnames(V)

    (result + t(result)) / 2
}



###### 13. Multinomial Probit Prediction Standard Errors ####

.econ_multinomial_probit_prediction_se = function(
    object,
    newdata,
    step = 1e-5
) {

    ###### 1. Construct prediction model matrix ####

    terms_object = stats::delete.response(object$terms)

    mf_new = stats::model.frame(
        terms_object,
        data = newdata,
        xlev = object$xlevels,
        na.action = stats::na.pass
    )

    X = stats::model.matrix(
        terms_object,
        data = mf_new,
        contrasts.arg = object$contrasts
    )

    C = stats::coef(object)

    X = X[, colnames(C), drop = FALSE]

    V = stats::vcov(object)

    beta = as.vector(t(C))

    parameter_names = unlist(
        lapply(rownames(C), function(category) {
            paste0(category, ":", colnames(C))
        }),
        use.names = FALSE
    )

    gamma = if (identical(object$error_covariance, "estimated")) {
        object$covariance_parameters
    } else {
        numeric(0)
    }
    theta = c(setNames(beta, parameter_names), gamma)

    if (!setequal(names(theta), colnames(V))) {
        stop(
            "Coefficient and covariance names differ.",
            call. = FALSE
        )
    }

    parameter_names = colnames(V)
    theta = theta[parameter_names]

    ###### 2. Prediction function ####

probabilities = function(par, X_new) {

    B = matrix(
        par[seq_along(beta)],
        nrow = nrow(C),
        ncol = ncol(C),
        byrow = TRUE
    )

    eta = cbind(
        0,
        X_new %*% t(B)
    )

    colnames(eta) = object$lev

    .econ_multinomial_probit_probs_general(
        eta = eta,
        normalize = TRUE,
        Omega = if (length(gamma) > 0L) {
            .econ_multinomial_probit_omega(
                par[names(gamma)],
                nrow(C)
            )
        } else {
            object$Omega
        }
    )
}

    ###### 3. Initialize standard errors ####

    result = matrix(
        NA_real_,
        nrow = nrow(X),
        ncol = length(object$lev),
        dimnames = list(NULL, object$lev)
    )

    complete = stats::complete.cases(X)

    if (!any(complete)) {
        return(result)
    }

    X_complete = X[complete, , drop = FALSE]

    n = nrow(X_complete)
    J = length(object$lev)
    K = length(theta)

    ###### 4. Numerical probability gradients ####

    gradients = array(
        0,
        dim = c(n, J, K)
    )

    for (k in seq_len(K)) {

        h = step * max(1, abs(theta[k]))

        theta_plus = theta
        theta_minus = theta

        theta_plus[k] = theta_plus[k] + h
        theta_minus[k] = theta_minus[k] - h

        P_plus = probabilities(theta_plus, X_complete)
        P_minus = probabilities(theta_minus, X_complete)

        gradients[, , k] = (
            P_plus - P_minus
        ) / (2 * h)
    }

    ###### 5. Delta-method variances ####

    for (j in seq_len(J)) {

        gradient = matrix(
            gradients[, j, ],
            nrow = n,
            ncol = K
        )

        variance = rowSums(
            (gradient %*% V) * gradient
        )

        result[complete, j] = sqrt(
            pmax(0, variance)
        )
    }

    result
}
