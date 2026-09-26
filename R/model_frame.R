.econ_model_frame <- function(
    formula,
    data,
    weights_expr = NULL,
    subset_expr = NULL,
    cluster_expr = NULL,
    na.action = stats::na.omit,
    enclos = parent.frame()
) {

    ###### 1. Validate data ####

    if (!is.data.frame(data)) {
        stop(
            "'data' must be a data.frame.",
            call. = FALSE
        )
    }

    #### Validate NA action ####

    if (is.character(na.action)) {

        if (length(na.action) != 1L) {
            stop(
                "'na.action' must specify a single NA-handling function.",
                call. = FALSE
            )
        }

        na.action <- get(
            na.action,
            mode = "function"
        )
    }

    if (!is.function(na.action)) {
        stop(
            "'na.action' must be a function or the name of a function.",
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


    ###### 3. Evaluate cluster inside data ####

    cluster_full <- NULL

    if (!is.null(cluster_expr)) {

        cluster_full <- eval(
            cluster_expr,
            envir = data,
            enclos = enclos
        )

        if (length(cluster_full) != nrow(data)) {
            stop(
                "'cluster' must have one value per row of data.",
                call. = FALSE
            )
        }
    }


    ###### 4. Evaluate subset inside data ####

    subset_value <- NULL

    if (!is.null(subset_expr)) {

        subset_value <- eval(
            subset_expr,
            envir = data,
            enclos = enclos
        )

        if (
            !is.logical(subset_value) &&
            !is.numeric(subset_value)
        ) {
            stop(
                "'subset' must evaluate to a logical or numeric vector.",
                call. = FALSE
            )
        }
    }


    ###### 5. Determine observations eligible for warnings ####

    mf_no_extra <- stats::model.frame(
        formula = formula,
        data = data,
        na.action = stats::na.pass,
        drop.unused.levels = TRUE
    )

    model_complete <- stats::complete.cases(
        mf_no_extra
    )

    subset_complete <- rep.int(
        TRUE,
        nrow(data)
    )

    if (!is.null(subset_value)) {

        if (is.numeric(subset_value)) {

            subset_complete <- seq_len(
                nrow(data)
            ) %in% subset_value

        } else {

            subset_complete <- !is.na(subset_value) &
                subset_value
        }
    }

    warning_sample <- model_complete &
        subset_complete


###### 6. Warn about missing weights ####

if (
    !identical(
        na.action,
        stats::na.fail
    ) &&
    !is.null(weights_expr) &&
    anyNA(w_full)
) {

    missing_weights_only <- is.na(w_full) &
        warning_sample

    if (any(missing_weights_only)) {
        warning(
            paste0(
                "Some observations have missingness in the weights variable(s) ",
                "but not in the outcome or covariates. ",
                "These observations have been dropped."
            ),
            call. = FALSE
        )
    }
}


###### 7. Warn about missing cluster IDs ####

if (
    !identical(
        na.action,
        stats::na.fail
    ) &&
    !is.null(cluster_expr) &&
    anyNA(cluster_full)
) {

    missing_cluster_only <- is.na(cluster_full) &
        warning_sample

    if (any(missing_cluster_only)) {
        warning(
            paste0(
                "Some observations have missingness in the cluster variable(s) ",
                "but not in the outcome or covariates. ",
                "These observations have been dropped."
            ),
            call. = FALSE
        )
    }
}


    ###### 8. Construct model frame ####

    args <- list(
        formula = formula,
        data = data,
        weights = w_full,
        na.action = stats::na.pass,
        drop.unused.levels = TRUE
    )

    if (!is.null(subset_value)) {
        args$subset <- subset_value
    }

    mf <- do.call(
        stats::model.frame,
        args
    )


    ###### 9. Determine original row positions ####

    used <- suppressWarnings(
        as.integer(row.names(mf))
    )

    if (anyNA(used)) {

        used <- match(
            row.names(mf),
            row.names(data)
        )
    }


    ###### 10. Add cluster to NA handling ####

    cluster <- NULL

    if (!is.null(cluster_full)) {

        cluster <- cluster_full[
            used
        ]

        complete <- stats::complete.cases(mf) &
            !is.na(cluster)

    } else {

        complete <- stats::complete.cases(mf)
    }


    ###### 11. Apply NA action ####

    if (!all(complete)) {

        if (identical(
            na.action,
            stats::na.fail
        )) {
            stop(
                "missing values in object",
                call. = FALSE
            )
        }

        mf <- mf[
            complete,
            ,
            drop = FALSE
        ]

        used <- used[
            complete
        ]

        if (!is.null(cluster)) {
            cluster <- cluster[
                complete
            ]
        }
    }


    ###### 12. Store omitted observations ####

    omitted <- which(
        !seq_len(nrow(data)) %in% used
    )

    if (length(omitted) == 0) {
        omitted <- NULL
    }


    ###### 13. Extract estimation weights ####

    w <- stats::model.weights(mf)

    if (is.null(w)) {
        w <- rep.int(
            1,
            nrow(mf)
        )
    }


    #### 14. Validate weights ####

    if (!is.numeric(w)) {
        stop(
            "'weights' must be a numeric vector.",
            call. = FALSE
        )
    }

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


    #### 15. Check for positive weights ####

    #if (!any(w > 0)) {
    #    stop(
    #        "At least one observation must have a positive weight.",
    #        call. = FALSE
    #    )
    #}


    #### 16. Return model-frame information ####

    list(
        frame = mf,
        used = used,
        omitted = omitted,
        weights = w,
        cluster = cluster
    )
}