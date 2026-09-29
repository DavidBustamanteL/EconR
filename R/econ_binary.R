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
#' @param link Binary-response link function. Currently "logit" is supported.
#'
#' @return An object of class "econ_binary".
#'
#' @export
econ_binary <- function(
    formula,
    data,
    weights = NULL,
    subset = NULL,
    na.action = stats::na.omit,
    link = c(
        "logit",
        "probit"
    )
) {

    ###### 1. Validate link ####

    link <- match.arg(
        link
    )

    if (link != "logit") {
        stop(
            "Only link = 'logit' is currently implemented.",
            call. = FALSE
        )
    }


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


    ###### 3. Construct estimation sample ####

    mf0 <- .econ_model_frame(
        formula = formula,
        data = data,
        weights_expr = weights_expr,
        subset_expr = subset_expr,
        cluster_expr = NULL,
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
        coefficients = fit$coefficients,
        fitted.values = fit$fitted.values,
        linear.predictors = fit$linear.predictors,
        residuals = fit$residuals,
        deviance = fit$deviance,
        null.deviance = fit$null.deviance,
        df.residual = fit$df.residual,
        rank = fit$rank,
        converged = fit$converged,
        iter = fit$iter,
        link = link,
        used = mf0$used,
        omitted = mf0$omitted
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