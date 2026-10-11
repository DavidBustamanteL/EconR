utils::globalVariables(c(
    ".econ_ordered_weights",
    ".econ_ordered_subset"
))

###### 1. Ordered Logit and Probit Estimation ####

#' Ordered Logit and Probit Regression
#'
#' @param formula Model formula.
#' @param data Data frame containing model variables.
#' @param link Either "logit" or "probit".
#' @param weights Optional observation weights.
#' @param subset Optional subset of observations.
#' @param na.action Missing-value handling function.
#' @param vcov Covariance estimator: "classical", "HC0", "HC1", or "CR0".
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
            stop("CR0 requires a cluster variable.")
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
