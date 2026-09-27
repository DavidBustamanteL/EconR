.econ_model_frame <- function(
    formula,
    data,
    weights_expr = NULL,
    subset_expr = NULL,
    cluster_expr = NULL,
    na.action = stats::na.omit,
    enclos = parent.frame()
) {

    #### 1. Validate data ####

    if (!is.data.frame(data)) {
        stop(
            "'data' must be a data.frame.",
            call. = FALSE
        )
    }


    #### 2. Evaluate weights inside data ####

    if (is.null(weights_expr)) {

        w_full <- rep.int(
            1,
            nrow(data)
        )

    } else {

        w_full <- eval(
            weights_expr,
            envir = data,
            enclos = enclos
        )

        if (length(w_full) != nrow(data)) {
            stop(
                "'weights' must have one value per row of data.",
                call. = FALSE
            )
        }
    }


    #### 3. Evaluate subset inside data ####

    subset_value <- NULL

    if (!is.null(subset_expr)) {

        subset_value <- eval(
            subset_expr,
            envir = data,
            enclos = enclos
        )
    }


    #### 4. Construct model frame ####

    args <- list(
        formula = formula,
        data = data,
        weights = w_full,
        na.action = na.action,
        drop.unused.levels = TRUE
    )

    if (!is.null(subset_value)) {
        args$subset <- subset_value
    }

    mf <- do.call(
        stats::model.frame,
        args
    )


    #### 5. Store omitted observations ####

    omitted <- attr(
        mf,
        "na.action"
    )


    #### 6. Determine observations used ####

    used <- suppressWarnings(
        as.integer(
            row.names(mf)
        )
    )

    if (anyNA(used)) {

        used <- match(
            row.names(mf),
            row.names(data)
        )
    }


    #### 7. Extract estimation weights ####

    w <- stats::model.weights(
        mf
    )

    if (is.null(w)) {

        w <- rep.int(
            1,
            nrow(mf)
        )
    }


    #### 8. Validate weights ####

    if (anyNA(w)) {
        stop(
            "Missing values remain in weights after NA handling.",
            call. = FALSE
        )
    }

    if (any(!is.finite(w))) {
        stop(
            "Weights must be finite.",
            call. = FALSE
        )
    }

    if (any(w < 0)) {
        stop(
            "Negative weights are not allowed.",
            call. = FALSE
        )
    }


    #### 9. Return model-frame information ####

    list(
        frame = mf,
        used = used,
        omitted = omitted,
        weights = w
    )
}


#' Fit an OLS or WLS model
#'
#' Uses standard formula/model-frame semantics: only observations missing values
#' required by the current specification are omitted. Robust covariance matrices
#' are provided by the sandwich and clubSandwich packages.
#'
#' @param formula model formula
#' @param data data.frame
#' @param weights optional non-negative observation weights
#' @param subset optional subset
#' @param na.action missing-data action; defaults to na.omit
#' @param vcov covariance type: classical, HC0, HC1, HC2, HC3, CR0, CR2, or stata
#' @param cluster optional cluster variable/vector; required for CR0, CR2, and stata
#'
#' @export
econ_lm <- function(
    formula,
    data,
    weights = NULL,
    subset = NULL,
    na.action = stats::na.omit,
    vcov = c(
        "classical",
        "HC0",
        "HC1",
        "HC2",
        "HC3",
        "CR0",
        "CR2",
        "stata"
    ),
    cluster = NULL
) {

    #### 1. Match covariance estimator ####

    vcov <- match.arg(
        vcov
    )


    #### 2. Capture original expressions ####

    weights_expr <- substitute(
        weights
    )

    subset_expr <- substitute(
        subset
    )

    cluster_expr <- substitute(
        cluster
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

    if (identical(
        cluster_expr,
        quote(NULL)
    )) {
        cluster_expr <- NULL
    }


    #### 3. Construct estimation sample ####

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

    if (nrow(mf) == 0) {
    stop(
        "No complete observations available for estimation.",
        call. = FALSE
    )
}
    
    if (!any(mf0$weights > 0)) {
    stop(
        "At least one observation must have a positive weight.",
        call. = FALSE
    )
}


    #### 4. Extract model components ####

    mt <- attr(
        mf,
        "terms"
    )

    y <- stats::model.response(
        mf
    )

    if (is.null(y)) {
        stop(
            "A response variable must be specified in the formula.",
            call. = FALSE
        )
    }

    if (is.matrix(y) || NCOL(y) != 1L) {
        stop(
            "Multivariate responses are not supported by econ_lm().",
            call. = FALSE
        )
    }

    if (!is.numeric(y)) {
        stop(
            "The response variable must be numeric.",
            call. = FALSE
        )
    }

    offset <- stats::model.offset(
        mf
    )

    if (is.null(offset)) {
        offset <- rep.int(
            0,
            nrow(mf)
        )
    }

    y_adjusted <- y - offset

    X <- stats::model.matrix(
        mt,
        mf
    )

    w <- mf0$weights

    #### 5. Construct weighted design ####

    sw <- sqrt(
        w
    )

    Xw <- X * sw

    yw <- y_adjusted * sw


    #### 6. Estimate coefficients ####

    fit <- stats::lm.fit(
        Xw,
        yw
    )

    beta <- fit$coefficients

    rank <- fit$rank

    estimable <- fit$qr$pivot[
        seq_len(rank)
    ]

    aliased <- setdiff(
        seq_len(ncol(X)),
        estimable
    )


    #### 7. Calculate fitted values and residuals ####

    fitted <- as.vector(
        X[, estimable, drop = FALSE] %*%
            beta[estimable] +
            offset
    )

    resid <- as.vector(
        y - fitted
    )

    names(fitted) <- row.names(mf)
    names(resid) <- row.names(mf)

    rank <- fit$rank

    df.residual <- sum(w > 0) - rank


#### 8. Check for essentially perfect fit ####

if (df.residual > 0) {

    rss <- sum(
        w * resid^2
    )

    y_center <- stats::weighted.mean(
        y,
        w
    )

    tss <- sum(
        w * (y - y_center)^2
    )

    scale <- max(
        sum(w * fitted^2),
        tss,
        .Machine$double.eps
    )

    if (rss < 1e-20 * scale) {
        warning(
            "Essentially perfect fit: inference may be unreliable.",
            call. = FALSE
        )
    }
}


    #### 9. Extract cluster variable ####

    cl_used <- mf0$cluster


    #### 10. Validate covariance specification ####

    cluster_types = c(
        "CR0",
        "CR2",
        "stata"
    )

    noncluster_types = c(
        "classical",
        "HC0",
        "HC1",
        "HC2",
        "HC3"
    )

    if (
        vcov %in% cluster_types &&
        is.null(cl_used)
    ) {
        stop(
            paste0(
                "vcov = '",
                vcov,
                "' requires a cluster variable."
            ),
            call. = FALSE
        )
    }

    if (
        !is.null(cl_used) &&
        vcov %in% noncluster_types
    ) {
        stop(
            paste0(
                "A cluster variable was supplied, but vcov = '",
                vcov,
                "' is not a clustered covariance estimator."
            ),
            call. = FALSE
        )
    }


    #### 11. Construct model object ####

    out <- list(
        call = match.call(),
        formula = formula,
        terms = mt,
        model = mf,
        x = X,
        y = y,
        weights = w,
        offset = offset,
        coefficients = beta,
        fitted.values = fitted,
        residuals = resid,
        rank = rank,
        estimable = estimable,
        aliased = aliased,
        df.residual = df.residual,
        nobs = nrow(X),
        used = mf0$used,
        na.action = if (
            identical(
                na.action,
                stats::na.exclude
            ) &&
            !is.null(mf0$omitted)
        ) {
            structure(
                mf0$omitted,
                names = row.names(data)[mf0$omitted],
                class = "exclude"
            )
        } else {
            mf0$omitted
        },
        na.action_type = na.action,
        vcov_type = vcov,
        cluster = cl_used,
        contrasts = attr(
            X,
            "contrasts"
        ),
        xlevels = stats::.getXlevels(
            mt,
            mf
        )
    )

    class(out) <- "econ_lm"


    #### 12. Calculate covariance matrix ####

    out$vcov <- econ_vcov(
        out,
        type = vcov
    )


    #### 13. Return model ####

    out
}


#' Covariance matrix for econ_lm
#'
#' @param object econ_lm object
#' @param type covariance type
#'
#' @export
econ_vcov <- function(
    object,
    type = object$vcov_type
) {

    #### 1. Identify estimable coefficients ####

    estimable <- object$estimable

    if (is.null(estimable)) {
        estimable <- seq_len(
            ncol(object$x)
        )
    }

    coefficient_names <- names(
        object$coefficients
    )

    expand_vcov <- function(V) {

        V_full <- matrix(
            NA_real_,
            nrow = length(coefficient_names),
            ncol = length(coefficient_names),
            dimnames = list(
                coefficient_names,
                coefficient_names
            )
        )

        V_full[estimable, estimable] <- as.matrix(V)

        V_full
    }


    #### 2. Construct estimable model object ####

    reduced <- object

    reduced$x <- object$x[
        ,
        estimable,
        drop = FALSE
    ]

    reduced$coefficients <- object$coefficients[
        estimable
    ]

    reduced$estimable <- seq_len(
        length(estimable)
    )

    reduced$aliased <- integer(0)


#### 3. Classical covariance ####

if (type == "classical") {

    Xw <- reduced$x * sqrt(
        reduced$weights
    )

    s2 <- sum(
        reduced$weights *
            reduced$residuals^2
    ) / reduced$df.residual

    qr_Xw <- qr(
        Xw
    )

    R <- qr.R(
        qr_Xw
    )

    R <- R[
        seq_len(reduced$rank),
        seq_len(reduced$rank),
        drop = FALSE
    ]

    V <- s2 * chol2inv(
        R
    )

    return(
        expand_vcov(V)
    )
}


    ###### 4. Heteroskedasticity-robust covariance ####

    if (type %in% c(
        "HC0",
        "HC1",
        "HC2",
        "HC3"
    )) {

        V <- sandwich::vcovHC(
            reduced,
            type = type
        )

        return(
            expand_vcov(V)
        )
    }


    ###### 5. CR0 clustered covariance ####

    if (type == "CR0") {

        if (is.null(object$cluster)) {
            stop(
                "CR0 clustered covariance requires a cluster variable.",
                call. = FALSE
            )
        }

        V <- clubSandwich::vcovCR(
            reduced,
            cluster = reduced$cluster,
            type = "CR0"
        )

        return(
            expand_vcov(V)
        )
    }


    ###### 6. CR2 clustered covariance ####

    if (type == "CR2") {

        if (is.null(object$cluster)) {
            stop(
                "CR2 clustered covariance requires a cluster variable.",
                call. = FALSE
            )
        }

        V <- clubSandwich::vcovCR(
            reduced,
            cluster = reduced$cluster,
            type = "CR2"
        )

        return(
            expand_vcov(V)
        )
    }


    ###### 7. Stata-style clustered covariance ####

    if (type == "stata") {

        if (is.null(object$cluster)) {
            stop(
                "Stata-style clustered covariance requires a cluster variable.",
                call. = FALSE
            )
        }

        V <- sandwich::vcovCL(
            reduced,
            cluster = reduced$cluster,
            type = "HC1"
        )

        return(
            expand_vcov(V)
        )
    }


    ###### 8. Unknown covariance estimator ####

    stop(
        "Unknown covariance type.",
        call. = FALSE
    )
}



###### 1. Sandwich compatibility ####

#' @export
bread.econ_lm <- function(
    x,
    ...
) {

    estimable <- x$estimable

    if (is.null(estimable)) {
        estimable <- seq_len(
            ncol(x$x)
        )
    }

    X <- x$x[
        ,
        estimable,
        drop = FALSE
    ]

    w <- x$weights
    n <- stats::nobs(x)

    n * solve(
        crossprod(
            X,
            X * w
        )
    )
}


#' @export
estfun.econ_lm <- function(
    x,
    ...
) {

    estimable <- x$estimable

    if (is.null(estimable)) {
        estimable <- seq_len(
            ncol(x$x)
        )
    }

    X <- x$x[
        ,
        estimable,
        drop = FALSE
    ]

    as.vector(
        x$weights *
            x$residuals
    ) * X
}


#' @export
hatvalues.econ_lm <- function(
    model,
    ...
) {

    estimable <- model$estimable

    if (is.null(estimable)) {
        estimable <- seq_len(
            ncol(model$x)
        )
    }

    Xw <- model$x[
        ,
        estimable,
        drop = FALSE
    ] * sqrt(
        model$weights
    )

    XtX_inv <- chol2inv(
        qr.R(
            qr(Xw)
        )
    )

    h <- rowSums(
        (Xw %*% XtX_inv) *
            Xw
    )

    names(h) <- rownames(
        model$x
    )

    h
}


#### 2. Standard model methods ####

#' @export
model.matrix.econ_lm <- function(
    object,
    ...
) {

    object$x
}


#' @export
weights.econ_lm <- function(
    object,
    ...
) {

    object$weights
}


#' @export
nobs.econ_lm <- function(
    object,
    ...
) {

    w <- stats::weights(object)

    if (is.null(w)) {
        return(
            object$nobs
        )
    }

    sum(w != 0)
}


#' @export
coef.econ_lm <- function(
    object,
    ...
) {

    object$coefficients
}


#' @export
vcov.econ_lm <- function(
    object,
    ...
) {

    object$vcov
}

#' @export
confint.econ_lm <- function(
    object,
    parm,
    level = 0.95,
    ...
) {

    ###### 1. Extract coefficients ####

    cf <- stats::coef(
        object
    )

    pnames <- names(
        cf
    )


    ###### 2. Select parameters ####

    if (missing(parm)) {

        parm <- pnames

    } else if (is.numeric(parm)) {

        parm <- pnames[
            parm
        ]
    }


    ###### 3. Compute confidence level ####

    alpha <- (
        1 - level
    ) / 2

    probs <- c(
        alpha,
        1 - alpha
    )

    critical <- stats::qt(
        probs,
        df = object$df.residual
    )


    ###### 4. Compute confidence intervals ####

    se <- sqrt(
        diag(
            stats::vcov(
                object
            )
        )
    )

    ci <- array(
        NA_real_,
        dim = c(
            length(parm),
            2L
        ),
        dimnames = list(
            parm,
            paste0(
                format(
                    100 * probs,
                    trim = TRUE
                ),
                " %"
            )
        )
    )

    ci[] <- cf[parm] +
        se[parm] %o%
        critical


    ###### 5. Return ####

    ci
}

#' @export
anova.econ_lm <- function(
    object,
    ...,
    scale = 0,
    test = "F"
) {

    ###### 1. Collect models ####

    objects <- list(
        object,
        ...
    )

    nmodels <- length(
        objects
    )


    ###### 2. Handle single model ####

    if (nmodels == 1L) {

        X <- object$x

        y <- object$y - object$offset

        w <- object$weights

        assign <- attr(
            X,
            "assign"
        )

        term_labels <- attr(
            object$terms,
            "term.labels"
        )


        ###### 3. Construct weighted estimation problem ####

        sqrt_w <- sqrt(
            w
        )

        Xw <- X * sqrt_w

        yw <- y * sqrt_w


        ###### 4. Compute sequential decomposition ####

        fit <- stats::lm.fit(
            Xw,
            yw
        )

        rank <- fit$rank

        pivot <- fit$qr$pivot[
            seq_len(rank)
        ]

        effects <- fit$effects[
            seq_len(rank)
        ]

        assign_estimable <- assign[
            pivot
        ]

        split_effects <- split(
            effects^2,
            assign_estimable
        )

        ss <- vapply(
            split_effects,
            sum,
            numeric(1)
        )

        df <- lengths(
            split(
                assign_estimable,
                assign_estimable
            )
        )


        ###### 5. Remove intercept ####

        term_ids <- as.integer(
            names(ss)
        )

        keep <- term_ids != 0L

        ss <- ss[
            keep
        ]

        df <- df[
            keep
        ]

        term_ids <- term_ids[
            keep
        ]

        labels <- term_labels[
            term_ids
        ]


        ###### 6. Compute residual quantities ####

        ssr <- sum(
            w *
                object$residuals^2
        )

        dfr <- object$df.residual

        residual_ms <- ssr / dfr

        ms <- ss / df

        f <- ms / residual_ms

        p <- stats::pf(
            f,
            df,
            dfr,
            lower.tail = FALSE
        )


        ###### 7. Construct single-model table ####

        table <- data.frame(
            Df = c(
                df,
                dfr
            ),
            `Sum Sq` = c(
                ss,
                ssr
            ),
            `Mean Sq` = c(
                ms,
                residual_ms
            ),
            `F value` = c(
                f,
                NA_real_
            ),
            `Pr(>F)` = c(
                p,
                NA_real_
            ),
            check.names = FALSE
        )

        rownames(
            table
        ) <- c(
            labels,
            "Residuals"
        )


        ###### 8. Return single-model ANOVA ####

        return(
            structure(
                table,
                heading = c(
                    "Analysis of Variance Table\n",
                    paste(
                        "Response:",
                        deparse(
                            object$formula[[2L]]
                        )
                    )
                ),
                class = c(
                    "anova",
                    "data.frame"
                )
            )
        )
    }


    ###### 9. Validate model classes ####

    valid <- vapply(
        objects,
        inherits,
        logical(1),
        what = "econ_lm"
    )

    if (!all(valid)) {
        stop(
            "All models must be 'econ_lm' objects.",
            call. = FALSE
        )
    }


    ###### 10. Check response variables ####

    responses <- vapply(
        objects,
        function(x) {
            paste(
                deparse(
                    x$terms[[2L]]
                ),
                collapse = ""
            )
        },
        character(1)
    )

    same_response <- responses == responses[1L]

    if (!all(same_response)) {

        warning(
            sprintf(
                "models with response %s removed because response differs from model 1",
                paste(
                    sQuote(
                        responses[
                            !same_response
                        ]
                    ),
                    collapse = ", "
                )
            ),
            call. = FALSE
        )

        objects <- objects[
            same_response
        ]
    }

    nmodels <- length(
        objects
    )

    if (nmodels == 1L) {
        return(
            anova(
                objects[[1L]]
            )
        )
    }


    ###### 11. Check estimation sample sizes ####

    ns <- vapply(
        objects,
        function(x) {
            length(
                x$residuals
            )
        },
        integer(1)
    )

    if (any(ns != ns[1L])) {
        stop(
            "models were not all fitted to the same size of dataset",
            call. = FALSE
        )
    }


    ###### 12. Extract residual quantities ####

    resdf <- vapply(
        objects,
        function(x) {
            x$df.residual
        },
        numeric(1)
    )

    resdev <- vapply(
        objects,
        function(x) {
            sum(
                x$weights *
                    x$residuals^2
            )
        },
        numeric(1)
    )


    ###### 13. Construct comparison table ####

    table <- data.frame(
        Res.Df = resdf,
        RSS = resdev,
        Df = c(
            NA_real_,
            -diff(
                resdf
            )
        ),
        `Sum of Sq` = c(
            NA_real_,
            -diff(
                resdev
            )
        ),
        check.names = FALSE
    )

    rownames(
        table
    ) <- as.character(
        seq_len(
            nmodels
        )
    )


    ###### 14. Compute F tests ####

    if (!is.null(test)) {

        bigmodel <- order(
            resdf
        )[1L]

        if (scale > 0) {

            scale_value <- scale

        } else {

            scale_value <- resdev[
                bigmodel
            ] / resdf[
                bigmodel
            ]
        }

        df_difference <- table$Df

        ss_difference <- table$`Sum of Sq`

        f <- (
            ss_difference /
                df_difference
        ) / scale_value

        p <- stats::pf(
            f,
            abs(
                df_difference
            ),
            resdf[
                bigmodel
            ],
            lower.tail = FALSE
        )

        table$F <- f

        table$`Pr(>F)` <- p
    }


    ###### 15. Construct comparison heading ####

    variables <- vapply(
        objects,
        function(x) {
            paste(
                deparse(
                    stats::formula(
                        x
                    )
                ),
                collapse = "\n"
            )
        },
        character(1)
    )

    topnote <- paste0(
        "Model ",
        format(
            seq_len(
                nmodels
            )
        ),
        ": ",
        variables,
        collapse = "\n"
    )


    ###### 16. Return model-comparison ANOVA ####

    structure(
        table,
        heading = c(
            "Analysis of Variance Table\n",
            topnote
        ),
        class = c(
            "anova",
            "data.frame"
        )
    )
}


#' @export
residuals.econ_lm <- function(
    object,
    ...
) {

    resid <- object$residuals

    if (
        identical(
            object$na.action_type,
            stats::na.exclude
        ) &&
        !is.null(object$na.action)
    ) {

        resid <- stats::napredict(
            object$na.action,
            resid
        )
    }

    resid
}


#' @export
fitted.econ_lm <- function(
    object,
    ...
) {

    fitted <- object$fitted.values

    if (
        identical(
            object$na.action_type,
            stats::na.exclude
        ) &&
        !is.null(object$na.action)
    ) {

        fitted <- stats::napredict(
            object$na.action,
            fitted
        )
    }

    fitted
}


###### 3. Prediction ####

#' @export
predict.econ_lm <- function(
    object,
    newdata = NULL,
    na.action = stats::na.pass,
    ...
) {

    if (is.null(newdata)) {
        return(
            stats::fitted(object)
        )
    }

    tt <- stats::delete.response(
        object$terms
    )

    mf <- stats::model.frame(
        tt,
        newdata,
        na.action = na.action,
        xlev = object$xlevels
    )

    offset <- stats::model.offset(
        mf
    )

    if (is.null(offset)) {
        offset <- rep.int(
            0,
            nrow(mf)
        )
    }

    X <- stats::model.matrix(
        tt,
        mf,
        contrasts.arg = object$contrasts
    )

    estimable <- object$estimable

    if (is.null(estimable)) {
        estimable <- seq_len(
            ncol(X)
        )
    }

    pred <- as.vector(
        X[, estimable, drop = FALSE] %*%
            object$coefficients[estimable] +
            offset
    )

    names(pred) <- row.names(mf)

    pred
}


###### 4. Summary ####

#' @export
summary.econ_lm <- function(
    object,
    ...
) {

    ###### 1. Extract estimates ####

    est <- object$coefficients


    ###### 2. Compute inference statistics ####

    if (object$df.residual > 0) {

        se <- sqrt(
            diag(
                object$vcov
            )
        )

        stat <- est / se

        p <- 2 * stats::pt(
            abs(stat),
            df = object$df.residual,
            lower.tail = FALSE
        )

    } else {

        se <- rep(
            NaN,
            length(est)
        )

        stat <- rep(
            NaN,
            length(est)
        )

        p <- rep(
            NaN,
            length(est)
        )

        names(se) <- names(est)
        names(stat) <- names(est)
        names(p) <- names(est)
    }


    ###### 3. Construct coefficient table ####

    z <- cbind(
        Estimate = est,
        `Std. Error` = se,
        `t value` = stat,
        `Pr(>|t|)` = p
    )


    ###### 4. Compute model-fit statistics ####

    w <- stats::weights(
        object
    )

    if (is.null(w)) {
        w <- rep.int(
            1,
            length(object$residuals)
        )
    }

    rss <- sum(
        w * object$residuals^2
    )

    residual_se <- if (object$df.residual > 0) {

        sqrt(
            rss / object$df.residual
        )

    } else {

        NaN
    }

    y_adjusted <- object$y -
        object$offset

    n <- stats::nobs(
        object
    )

    if (attr(object$terms, "intercept") == 1L) {

        y_mean <- stats::weighted.mean(
            y_adjusted,
            w
        )

        tss <- sum(
            w *
                (y_adjusted - y_mean)^2
        )

    } else {

        tss <- sum(
            w *
                y_adjusted^2
        )
    }

    r_squared <- 1 -
        rss / tss

    if (object$df.residual > 0) {

        if (attr(object$terms, "intercept") == 1L) {

            adj_r_squared <- 1 -
                (1 - r_squared) *
                (
                    (n - 1) /
                        object$df.residual
                )

        } else {

            adj_r_squared <- 1 -
                (1 - r_squared) *
                (
                    n /
                        object$df.residual
                )
        }

    } else {

        adj_r_squared <- NaN
    }

    ###### 5. Compute joint model test ####

    coef_names <- names(
        est
    )

    test_terms <- which(
        coef_names != "(Intercept)" &
            !is.na(est)
    )

    if (
        length(test_terms) > 0L &&
        object$df.residual > 0
    ) {

        beta_test <- est[
            test_terms
        ]

        V_test <- object$vcov[
            test_terms,
            test_terms,
            drop = FALSE
        ]

        q <- length(
            beta_test
        )

        wald <- as.numeric(
            t(beta_test) %*%
                solve(
                    V_test,
                    beta_test
                )
        )

        f_statistic <- wald / q

        f_df1 <- q
        f_df2 <- object$df.residual

        f_p_value <- stats::pf(
            f_statistic,
            df1 = f_df1,
            df2 = f_df2,
            lower.tail = FALSE
        )

    } else {

        f_statistic <- NaN
        f_df1 <- 0L
        f_df2 <- object$df.residual
        f_p_value <- NaN
    }


    ###### 6. Construct summary object ####

    ans <- list(
        call = object$call,
        coefficients = z,
        df.residual = object$df.residual,
        nobs = object$nobs,
        vcov_type = object$vcov_type,
        sigma = residual_se,
        r.squared = r_squared,
        adj.r.squared = adj_r_squared,
        fstatistic = f_statistic,
        f.df1 = f_df1,
        f.df2 = f_df2,
        f.p.value = f_p_value
    )

    class(ans) <- "summary.econ_lm"

    ans
}


###### 5. Print summary ####

#' @export
print.summary.econ_lm <- function(
    x,
    ...
) {

    ###### 1. Print call ####

    cat(
        "\nCall:\n"
    )

    print(
        x$call
    )


    ###### 2. Print covariance type ####

    cat(
        "\nCovariance type: ",
        x$vcov_type,
        "\n",
        sep = ""
    )


    ###### 3. Print coefficients ####

    cat(
        "\nCoefficients:\n"
    )

    coef_table <- x$coefficients

    formatted <- matrix(
        "",
        nrow = nrow(coef_table),
        ncol = ncol(coef_table),
        dimnames = dimnames(coef_table)
    )

    formatted[, "Estimate"] <- sprintf(
        "%.5f",
        coef_table[, "Estimate"]
    )

    formatted[, "Std. Error"] <- sprintf(
        "%.5f",
        coef_table[, "Std. Error"]
    )

    formatted[, "t value"] <- sprintf(
        "%.5f",
        coef_table[, "t value"]
    )

    formatted[, "Pr(>|t|)"] <- sprintf(
        "%.5f",
        coef_table[, "Pr(>|t|)"]
    )

    print(
        formatted,
        quote = FALSE,
        right = TRUE
    )


    ###### 4. Print model-fit statistics ####

    cat(
        "\nResidual standard error: ",
        sprintf(
            "%.5f",
            x$sigma
        ),
        " on ",
        x$df.residual,
        " degrees of freedom",
        "\n",
        sep = ""
    )

    cat(
        "Multiple R-squared:  ",
        sprintf(
            "%.5f",
            x$r.squared
        ),
        "\n",
        sep = ""
    )

    cat(
        "Adjusted R-squared:  ",
        sprintf(
            "%.5f",
            x$adj.r.squared
        ),
        "\n",
        sep = ""
    )


    ###### 5. Print joint model test ####

    cat(
        "F-statistic: ",
        sprintf(
            "%.5f",
            x$fstatistic
        ),
        " on ",
        x$f.df1,
        " and ",
        x$f.df2,
        " DF, p-value: ",
        sprintf(
            "%.5f",
            x$f.p.value
        ),
        "\n",
        sep = ""
    )


    ###### 6. Print observations ####

    cat(
        "Observations: ",
        x$nobs,
        "\n",
        sep = ""
    )


    ###### 7. Return ####

    invisible(x)
}


###### 6. Print model ####

#' @export
print.econ_lm <- function(
    x,
    ...
) {

    print(
        summary(x)
    )

    invisible(x)
}


###### 7. Package Compatibility ####

#' Tidy an econ_lm model
#'
#' @param x econ_lm object
#' @param conf.int logical; include confidence intervals
#' @param conf.level confidence level
#' @param ... additional arguments
#'
#' @export
tidy.econ_lm <- function(
    x,
    conf.int = TRUE,
    conf.level = 0.95,
    ...
) {

    ###### 1. Extract coefficient table ####

    s <- summary(x)

    out <- data.frame(
        term = rownames(s$coefficients),
        estimate = s$coefficients[, "Estimate"],
        std.error = s$coefficients[, "Std. Error"],
        statistic = s$coefficients[, "t value"],
        p.value = s$coefficients[, "Pr(>|t|)"],
        row.names = NULL
    )


    ###### 2. Add confidence intervals ####

    if (conf.int) {

        if (x$df.residual > 0) {

            critical <- stats::qt(
                1 - (1 - conf.level) / 2,
                df = x$df.residual
            )

            out$conf.low <- out$estimate - critical * out$std.error
            out$conf.high <- out$estimate + critical * out$std.error

        } else {

            out$conf.low <- NaN
            out$conf.high <- NaN
        }
    }


    ###### 3. Add model information ####

    out$df <- x$df.residual

    out$outcome <- as.character(
        formula(x)[[2]]
    )


    ###### 4. Return tidy output ####

    out
}


#' Glance at an econ_lm model
#'
#' @param x econ_lm object
#' @param ... additional arguments
#'
#' @export
glance.econ_lm <- function(
    x,
    ...
) {

    ###### 1. Extract model information ####

    y <- stats::model.response(
        x$model
    )

    resid <- stats::residuals(
        x
    )

    w <- stats::weights(
        x
    )

    n <- stats::nobs(
        x
    )

    p <- x$rank


    ###### 2. Goodness-of-fit statistics ####

    rss <- sum(
        w * resid^2
    )

    has_intercept <- attr(
        stats::terms(x$model),
        "intercept"
    ) == 1

    if (has_intercept) {

        y_center <- stats::weighted.mean(
            y,
            w
        )

        tss <- sum(
            w * (y - y_center)^2
        )

        r_squared <- 1 - rss / tss

        if (x$df.residual > 0) {

            adj_r_squared <- 1 -
                (1 - r_squared) *
                (n - 1) /
                x$df.residual

        } else {

            adj_r_squared <- NaN
        }

    } else {

        tss <- sum(
            w * y^2
        )

        r_squared <- 1 - rss / tss

        if (x$df.residual > 0) {

            adj_r_squared <- 1 -
                (1 - r_squared) *
                n /
                x$df.residual

        } else {

            adj_r_squared <- NaN
        }
    }

    rmse <- sqrt(
        rss / sum(w)
    )


    ###### 3. Information criteria ####

    positive_weights <- w > 0

    log_lik <- -0.5 * (
        n * (
            log(2 * pi) +
            1 +
            log(rss / n)
        ) -
        sum(
            log(w[positive_weights])
        )
    )

    k <- p + 1

    aic <- -2 * log_lik +
        2 * k

    bic <- -2 * log_lik +
        log(n) * k


    ###### 4. Return glance output ####

    data.frame(
        r.squared = r_squared,
        adj.r.squared = adj_r_squared,
        aic = aic,
        bic = bic,
        rmse = rmse,
        nobs = n
    )
}