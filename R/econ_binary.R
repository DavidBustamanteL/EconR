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





utils::globalVariables(".econ_prior_weights")