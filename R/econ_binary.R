###### 1. Binary Response Estimator ####

#' Estimate a Binary Response Model
#'
#' Fits binary-response models using logit or probit links.
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param weights Optional observation weights.
#' @param subset Optional subset of observations.
#' @param na.action Function determining how missing values are handled.
#' @param link Binary-response link function. Either "logit" or "probit".
#' @param vcov Covariance estimator. Available options are "classical",
#'   "HC0", "HC1", "HC2", "HC3", "CR0", "CR2", and "stata".
#' @param cluster Optional clustering variable, evaluated within data.
#'   Required when vcov is "CR0", "CR2", or "stata".
#'
#' @return An object of class "econ_binary".
#'
#' @export
#' 
#' 
econ_binary <- function(
    formula,
    data,
    weights = NULL,
    subset = NULL,
    na.action = stats::na.omit,
    link = c("logit", "probit"),
    vcov = c(
        "classical", "HC0", "HC1", "HC2", "HC3",
        "CR0", "CR2", "stata"
    ),
    cluster = NULL
) {

    ###### 1. Validate link ####

    link = match.arg(link)

    vcov = match.arg(vcov)


    ###### 2. Capture expressions ####

    weights_expr <- substitute(
        weights
    )

    subset_expr <- substitute(
        subset
    )

    if (identical(
        weights_expr,
        quote(NULL)
    )) {
        weights_expr <- NULL
    }

    if (identical(
        subset_expr,
        quote(NULL)
    )) {
        subset_expr <- NULL
    }

    cluster_expr = substitute(cluster)

    if (identical(cluster_expr, quote(NULL))) {
        cluster_expr = NULL
    }


    ###### 3. Construct estimation sample ####

    mf0 <- .econ_model_frame(
        formula = formula,
        data = data,
        weights_expr = weights_expr,
        subset_expr = subset_expr,
        cluster_expr = cluster_expr,
        na.action = na.action,
        enclos = parent.frame()
    )

    mf <- mf0$frame

    if (nrow(mf) == 0L) {
        stop(
            "No complete observations available for estimation.",
            call. = FALSE
        )
    }


    ###### 4. Extract model components ####

    mt <- attr(
        mf,
        "terms"
    )

    y <- stats::model.response(
        mf
    )

    X <- stats::model.matrix(
        mt,
        mf
    )

    w <- mf0$weights


    ###### 5. Validate response ####

    if (is.matrix(y) || NCOL(y) != 1L) {
        stop(
            "Multivariate responses are not supported by econ_binary().",
            call. = FALSE
        )
    }

    if (!is.numeric(y)) {
        stop(
            "The response variable must be numeric.",
            call. = FALSE
        )
    }

    if (!all(y %in% c(0, 1))) {
        stop(
            "The response variable must contain only 0 and 1.",
            call. = FALSE
        )
    }


    ###### 6. Estimate model ####

    fit <- stats::glm.fit(
        x = X,
        y = y,
        weights = w,
        family = stats::binomial(
            link = link
        )
    )


    ###### 7. Construct model object ####

    ans <- list(
        call = match.call(),
        formula = formula,
        terms = mt,
        model = mf,
        x = X,
        y = y,
        weights = w,
        working.weights = fit$weights,
        vcov_type = vcov,
        cluster = mf0$cluster,
        coefficients = fit$coefficients,
        fitted.values = fit$fitted.values,
        linear.predictors = fit$linear.predictors,
        residuals = fit$residuals,
        deviance = fit$deviance,
        null.deviance = fit$null.deviance,
        df.residual = fit$df.residual,
        rank = fit$rank,
        qr = fit$qr,
        family = fit$family,
        converged = fit$converged,
        iter = fit$iter,
        link = link,
        used = mf0$used,
        omitted = mf0$omitted,
        contrasts = attr(X, "contrasts"),
        xlevels = stats::.getXlevels(mt, mf)
    )

    class(ans) <- "econ_binary"

    ans
}


###### 2. Basic Methods ####

#' @export
coef.econ_binary <- function(
    object,
    ...
) {

    object$coefficients
}


#' @export
print.econ_binary <- function(
    x,
    ...
) {

    cat(
        "\nCall:\n"
    )

    print(
        x$call
    )

    cat(
        "\nCoefficients:\n"
    )

    print(
        stats::coef(x)
    )

    invisible(x)
}


###### 3. Binary Covariance Method ####

#' @export
vcov.econ_binary <- function(
    object,
    type = object$vcov_type,
    ...
) {

    ###### 1. Validate covariance type ####

    type = match.arg(
        type,
        c(
            "classical", "HC0", "HC1", "HC2", "HC3",
            "CR0", "CR2", "stata"
        )
    )

    ###### 2. Heteroskedasticity-Robust Covariance ####

    if (type %in% c("HC0", "HC1", "HC2", "HC3")) {

        return(
            sandwich::vcovHC(
                object,
                type = type,
                ...
            )
        )
    }

    ###### 3. CR0 Clustered Covariance ####

    if (type == "CR0") {

        if (is.null(object$cluster)) {
            stop(
                "CR0 clustered covariance requires a cluster variable.",
                call. = FALSE
            )
        }

        return(
            sandwich::vcovCL(
                object,
                cluster = object$cluster,
                type = "HC",
                cadjust = FALSE
            )
        )
    }


###### 4. CR2 Clustered Covariance ####

if (type == "CR2") {

    if (is.null(object$cluster)) {
        stop(
            "CR2 clustered covariance requires a cluster variable.",
            call. = FALSE
        )
    }

    glm_data = object$model
    glm_data$.econ_prior_weights = object$weights

    glm_ref = stats::glm(
        formula = object$formula,
        data = glm_data,
        weights = .econ_prior_weights,
        family = object$family
    )

    return(
        as.matrix(
            clubSandwich::vcovCR(
                glm_ref,
                cluster = object$cluster,
                type = "CR2"
            )
        )
    )
}


        ###### 5. Stata-Style Clustered Covariance ####

    if (type == "stata") {

        if (is.null(object$cluster)) {
            stop(
                "Stata-style clustered covariance requires a cluster variable.",
                call. = FALSE
            )
        }

        return(
            sandwich::vcovCL(
                object,
                cluster = object$cluster,
                type = "HC1"
            )
        )
    }

    ###### 6. Extract QR decomposition ####

    p = object$rank

    R = qr.R(object$qr)
    R = R[seq_len(p), seq_len(p), drop = FALSE]

    ###### 7. Compute classical covariance ####

    V = chol2inv(R)

    ###### 8. Restore coefficient ordering ####

    pivot = object$qr$pivot[seq_len(p)]

    V_full = matrix(
        NA_real_,
        nrow = length(object$coefficients),
        ncol = length(object$coefficients),
        dimnames = list(
            names(object$coefficients),
            names(object$coefficients)
        )
    )

    V_full[pivot, pivot] = V

    V_full
}


###### 4. Binary Model Summary ####

#' @export
summary.econ_binary <- function(
    object,
    ...
) {

    ###### 1. Extract coefficients ####

    estimate = stats::coef(
        object
    )

    std_error = sqrt(
        diag(
            stats::vcov(object)
        )
    )

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
        deviance = object$deviance,
        null.deviance = object$null.deviance,
        df.residual = object$df.residual,
        converged = object$converged,
        iter = object$iter
    )

    class(ans) = "summary.econ_binary"

    ans
}


###### 5. Print Binary Model Summary ####

#' @export
print.summary.econ_binary <- function(
    x,
    ...
) {

    ###### 1. Print call ####

    cat("\nCall:\n")

    print(x$call)


    ###### 2. Print model information ####

    cat(
        "\nBinary response model: ",
        x$link,
        "\n",
        sep = ""
    )


    ###### 3. Print coefficients ####

    cat("\nCoefficients:\n")

    coef_table = x$coefficients

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

    print(
        formatted,
        quote = FALSE,
        right = TRUE
    )


    ###### 4. Print deviance ####

    cat(
        "\nNull deviance: ",
        sprintf("%.5f", x$null.deviance),
        "\n",
        sep = ""
    )

    cat(
        "Residual deviance: ",
        sprintf("%.5f", x$deviance),
        " on ",
        x$df.residual,
        " degrees of freedom\n",
        sep = ""
    )


    ###### 5. Print convergence ####

    cat(
        "Converged: ",
        x$converged,
        " (",
        x$iter,
        " iterations)\n",
        sep = ""
    )

    invisible(x)
}


###### 6. Log-Likelihood Method ####

#' @export
logLik.econ_binary <- function(
    object,
    ...
) {

    ###### 1. Extract model components ####

    y = object$y

    eta = object$linear.predictors

    w = object$weights


    ###### 2. Compute stable log probabilities ####

    log_p = stats::binomial(
        link = object$link
    )$linkinv(eta)

    p = log_p

    loglik = sum(
        w * stats::dbinom(
            y,
            size = 1,
            prob = p,
            log = TRUE
        )
    )


    ###### 3. Construct logLik object ####

    structure(
        loglik,
        class = "logLik",
        df = object$rank,
        nobs = sum(w > 0)
    )
}


###### 7. Prediction Method ####

#' @export
predict.econ_binary <- function(
    object,
    newdata = NULL,
    type = c("link", "response"),
    ...
) {

    ###### 1. Validate prediction type ####

    type = match.arg(type)


    ###### 2. In-sample prediction ####

    if (is.null(newdata)) {

        eta = object$linear.predictors

    } else {

        ###### 3. Construct new model matrix ####

        mt = stats::delete.response(
            object$terms
        )

        mf = stats::model.frame(
            mt,
            data = newdata,
            na.action = stats::na.pass,
            xlev = object$xlevels
        )

        X = stats::model.matrix(
            mt,
            mf,
            contrasts.arg = object$contrasts
        )

        ###### 4. Compute linear predictor ####

        eta = drop(
            X %*% object$coefficients
        )
    }


    ###### 5. Return requested prediction ####

    if (type == "link") {
        return(eta)
    }

    stats::binomial(
        link = object$link
    )$linkinv(eta)
}


###### 8. Estimating Functions ####

#' @export
estfun.econ_binary <- function(x, ...) {

    ###### 1. Extract model components ####

    X = x$x

    ###### 2. Compute weighted working residuals ####

    wres = unname(
        as.vector(x$residuals) *
            x$working.weights
    )

    ###### 3. Construct observation-level scores ####

    rval = wres * X

    ###### 4. Remove model-matrix attributes ####

    attr(rval, "assign") = NULL
    attr(rval, "contrasts") = NULL

    rval
}


###### 9. Bread Matrix ####

#' @export
bread.econ_binary <- function(x, ...) {

    ###### 1. Extract classical covariance ####

    V = stats::vcov(
        x,
        type = "classical"
    )

    ###### 2. Apply sample-size scaling ####

    n = sum(x$weights > 0)

    V * n
}


###### 10. Model Matrix Method ####

#' @export
model.matrix.econ_binary <- function(object, ...) {

    object$x
}


###### 11. Leverage Values ####

#' @export
hatvalues.econ_binary <- function(model, ...) {

    ###### 1. Extract weighted design matrix ####

    X = model$x
    w = model$working.weights

    ###### 2. Compute leverage from weighted QR ####

    Xw = X * sqrt(w)

    Q = qr.Q(
        qr(Xw),
        complete = FALSE
    )

    ###### 3. Return diagonal of hat matrix ####

    h = rowSums(Q[, seq_len(model$rank), drop = FALSE]^2)

    names(h) = rownames(X)

    h
}


###### Number of Observations ####

#' @export
nobs.econ_binary <- function(object, ...) {
    sum(object$weights > 0)
}


###### Broom Compatibility Methods ####

#' Tidy a Binary Response Model
#'
#' @param x An econ_binary object.
#' @param conf.int Include confidence intervals.
#' @param conf.level Confidence level.
#' @param ... Additional arguments.
#'
#' @return A data frame containing coefficient estimates and inference statistics.
#'
#' @export
tidy.econ_binary <- function(
    x,
    conf.int = TRUE,
    conf.level = 0.95,
    ...
) {

    ###### 1. Extract coefficient table ####

    s = summary(x)

    out = data.frame(
        term = rownames(s$coefficients),
        estimate = s$coefficients[, 1],
        std.error = s$coefficients[, 2],
        statistic = s$coefficients[, 3],
        p.value = s$coefficients[, 4],
        row.names = NULL
    )


    ###### 2. Add confidence intervals ####

    if (conf.int) {

        if (
            !is.numeric(conf.level) ||
            length(conf.level) != 1L ||
            is.na(conf.level) ||
            conf.level <= 0 ||
            conf.level >= 1
        ) {
            stop(
                "conf.level must be between 0 and 1.",
                call. = FALSE
            )
        }

        critical = stats::qnorm(
            1 - (1 - conf.level) / 2
        )

        out$conf.low = out$estimate - critical * out$std.error
        out$conf.high = out$estimate + critical * out$std.error
    }


    ###### 3. Add model information ####

    out$df = x$df.residual

    out$outcome = as.character(
        x$formula[[2]]
    )


    ###### 4. Return tidy output ####

    out
}

#' @export
glance.econ_binary <- function(
    x,
    ...
) {
    ###### 1. Extract model information ####

    n = stats::nobs(x)

    log_lik = as.numeric(
        stats::logLik(x)
    )


    ###### 2. Goodness-of-fit statistics ####

    deviance = x$deviance

    null_deviance = x$null.deviance

    pseudo_r_squared = if (null_deviance > 0) {
        1 - deviance / null_deviance
    } else {
        NaN
    }


    ###### 3. Information criteria ####

    aic = stats::AIC(x)

    bic = stats::BIC(x)


    ###### 4. Return glance output ####

    data.frame(
        pseudo.r.squared = pseudo_r_squared,
        logLik = log_lik,
        aic = aic,
        bic = bic,
        deviance = deviance,
        null.deviance = null_deviance,
        df.residual = x$df.residual,
        nobs = n,
        converged = x$converged
    )
}

###### Confidence Intervals ####

#' @export
confint.econ_binary = function(object, parm, level = 0.95, ...) {

    if (!is.numeric(level) || length(level) != 1L ||
        is.na(level) || level <= 0 || level >= 1) {
        stop("'level' must be a single number between 0 and 1.")
    }

    b = coef(object)
    se = sqrt(diag(vcov(object, ...)))

    if (missing(parm)) {
        parm = seq_along(b)
    }

    if (is.character(parm)) {
        parm = match(parm, names(b))

        if (anyNA(parm)) {
            stop("Unknown coefficient name in 'parm'.")
        }
    }

    if (anyNA(parm) || any(parm < 1 | parm > length(b))) {
        stop("Invalid coefficient selection in 'parm'.")
    }

    z = stats::qnorm((1 + level) / 2)

    ci = cbind(
        b[parm] - z * se[parm],
        b[parm] + z * se[parm]
    )

    colnames(ci) = c(
        paste0(format(100 * (1 - level) / 2), " %"),
        paste0(format(100 * (1 + level) / 2), " %")
    )

    ci
}


###### Fitted Values ####

#' @export
fitted.econ_binary = function(object, ...) {

    object$fitted.values
}




###### Residuals ####

#' @export
residuals.econ_binary = function(
    object,
    type = c("deviance", "pearson", "response", "working", "partial"),
    ...
) {

    type = match.arg(type)

    y = object$y
    mu = object$fitted.values
    w = object$weights

    residuals = switch(
        type,

        response = y - mu,

        pearson = (y - mu) * sqrt(w) /
            sqrt(object$family$variance(mu)),

        working = object$residuals,

        deviance = {
            dev = object$family$dev.resids(y, mu, w)
            sign(y - mu) * sqrt(pmax(dev, 0))
        },

        partial = {

            glm_data = object$model
            glm_data$.econ_prior_weights = object$weights

            glm_ref = stats::glm(
                formula = object$formula,
                data = glm_data,
                weights = .econ_prior_weights,
                family = object$family
            )

            stats::residuals(glm_ref, type = "partial")
        }
    )

    residuals
}



###### Model Weights ####

#' @export
weights.econ_binary = function(object, type = c("prior", "working"), ...) {

    type = match.arg(type)

    switch(
        type,
        prior = object$weights,
        working = object$working.weights
    )
}



###### Model Formula ####

#' @export
formula.econ_binary = function(x, ...) {

    x$formula
}



###### Analysis of Deviance ####

#' @export
anova.econ_binary = function(
    object,
    ...,
    test = "Chisq"
) {

    glm_data = object$model
    glm_data$.econ_prior_weights = object$weights

    glm_ref = stats::glm(
        formula = object$formula,
        data = glm_data,
        weights = .econ_prior_weights,
        family = object$family
    )

    stats::anova(
        glm_ref,
        ...,
        test = test
    )
}




###### Analysis of Deviance ####

#' @export
anova.econ_binary = function(object, ..., test = "Chisq") {

    models = list(object, ...)

    if (!all(vapply(models, inherits, logical(1), "econ_binary"))) {
        stop("All supplied models must be 'econ_binary' objects.")
    }

    if (length(models) > 1L) {

        reference = models[[1L]]

        for (i in seq_along(models)[-1L]) {

            current = models[[i]]

            # Check estimation samples
            if (!identical(reference$used, current$used)) {
                stop("Models must use identical estimation samples.")
            }

            # Check observation weights
            if (!isTRUE(all.equal(
                unname(reference$weights),
                unname(current$weights)
            ))) {
                stop("Models must use identical observation weights.")
            }

            # Check dependent variable
            if (!isTRUE(all.equal(
                unname(reference$y),
                unname(current$y)
            ))) {
                stop("Models must use identical dependent-variable values.")
            }

            # Check shared regressors
            shared = intersect(
                colnames(reference$x),
                colnames(current$x)
            )

            if (length(shared) > 0L) {

                if (!isTRUE(all.equal(
                    unname(reference$x[, shared, drop = FALSE]),
                    unname(current$x[, shared, drop = FALSE])
                ))) {
                    stop("Models contain different values for shared regressors.")
                }
            }

            # Check link function
            if (!identical(reference$link, current$link)) {
                stop("Models must use the same link function.")
            }
        }
    }

    # Reconstruct equivalent glm objects
    glm_models = lapply(models, function(x) {

        glm_data = x$model
        glm_data$.econ_prior_weights = x$weights

        stats::glm(
            formula = x$formula,
            data = glm_data,
            weights = .econ_prior_weights,
            family = x$family
        )
    })

    # Single-model analysis of deviance
    if (length(glm_models) == 1L) {
        return(stats::anova(glm_models[[1L]], test = test))
    }

    # Multiple-model comparison
    do.call(
        stats::anova,
        c(glm_models, list(test = test))
    )
}



###### Model Augmentation ####

#' @export
augment.econ_binary = function(
    x,
    data = NULL,
    newdata = NULL,
    se_fit = FALSE,
    ...
) {

    # Reconstruct the equivalent glm for standard diagnostics
    glm_data = x$model
    glm_data$.econ_prior_weights = x$weights

    glm_ref = stats::glm(
        formula = x$formula,
        data = glm_data,
        weights = .econ_prior_weights,
        family = x$family
    )

    # Obtain standard broom augmentation
    if (!is.null(newdata)) {

        result = broom::augment(
            glm_ref,
            newdata = newdata,
            se_fit = FALSE,
            ...
        )

        prediction_data = newdata

    } else {

        result = broom::augment(
            glm_ref,
            se_fit = FALSE,
            ...
        )

        prediction_data = NULL
    }

    
    # Remove the artificial weights column for unweighted models
    if ("(weights)" %in% names(result) &&
        all(x$weights == 1)) {

        result[["(weights)"]] = NULL
    }

    
    # Remove existing prediction columns before inserting EconR values
    result = result[
        ,
        !names(result) %in% c(".fitted", ".se.fit"),
        drop = FALSE
    ]

    # Replace fitted values with EconR predictions
    result$.fitted = unname(
        predict(x, newdata = prediction_data, type = "link")
    )

    # Standard errors based on EconR's selected covariance matrix
    if (isTRUE(se_fit)) {

        if (is.null(prediction_data)) {

            X = model.matrix(x)

        } else {

            terms_x = stats::delete.response(x$terms)

            X = stats::model.matrix(
                terms_x,
                data = prediction_data,
                contrasts.arg = x$contrasts,
                xlev = x$xlevels
            )
        }

        V = vcov(x)

        se = sqrt(pmax(
            rowSums((X %*% V) * X),
            0
        ))

        result$.se.fit = unname(se)

        

        # Place prediction columns before diagnostic columns
        prediction_columns = c(".fitted", ".se.fit")

        diagnostic_columns = c(
            ".resid",
            ".hat",
            ".sigma",
            ".cooksd",
            ".std.resid"
        )

        other_columns = setdiff(
            names(result),
            c(prediction_columns, diagnostic_columns)
        )

        column_order = c(
            other_columns,
            intersect(prediction_columns, names(result)),
            intersect(diagnostic_columns, names(result))
        )

        result = result[, column_order, drop = FALSE]

    } 

    result

} 










utils::globalVariables(".econ_prior_weights")