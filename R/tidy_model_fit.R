#' Extract variable importance (betas) from glmnet model
#' @param x A glmnet model object
#' @param phecode_description A logical indicating whether to include the phecode description
#' @return A data table with the exposure and beta
#' @importFrom data.table data.table
#' @importFrom dplyr filter mutate arrange desc left_join n select
#' @importFrom stringr str_wrap
#' @importFrom ms pheinfox
#' @importFrom stats coef
#' @export
tidy_glmnet_betas <- function(x, phecode_description = TRUE) {
    data.table::data.table(
        exposure = rownames(stats::coef(x)),
        beta = stats::coef(x)[,1]
    ) |>
    dplyr::filter(exposure != "(Intercept)") |>
    dplyr::mutate(
        abs_beta = abs(beta),
        rel_beta = abs_beta / max(abs_beta, na.rm = TRUE),
        sign = ifelse(beta < 0, "Negative", "Positive")
    ) |>
    dplyr::arrange(dplyr::desc(abs_beta)) |>
    dplyr::mutate(rank = 1:dplyr::n()) |>
    (\(x) if (phecode_description) {
        x |>
            dplyr::left_join(
                ms::pheinfox |> dplyr::select(exposure = phecode, description),
                by = "exposure"
            ) |>
            dplyr::mutate(
                label = stringr::str_wrap(paste0(description, " [", exposure, "]"), width = 35)
            )
    })()
}

#' Plot variable importance (betas) from glmnet model
#' @param x A glmnet model object
#' @param exposure_var The name of the exposure variable
#' @param beta_var The name of the beta variable to plot
#' @param top_n The number of top variables to plot
#' @return A ggplot2 object
#' @importFrom dplyr slice_head
#' @importFrom ggplot2 ggplot aes geom_bar labs coord_flip
#' @importFrom stats reorder
#' @export
plot_glmnet_vip <- function(x, exposure_var = label, beta_var = rel_beta, top_n = 15) {
    x |>
        dplyr::slice_head(n = top_n) |>
        ggplot2::ggplot(ggplot2::aes(x = stats::reorder(label, {{ beta_var }}), y = {{ beta_var }}, fill = sign)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
            x = "Variables",
            y = "Relative importance",
            title = paste0("Variable Importance (top ", top_n, ")")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::coord_flip()
}

#' Fit a glmnet model and extract variable importance (betas) and plot
#' @param data A data frame
#' @param exposures A character vector of exposure variables
#' @param outcome A character vector of outcome variables
#' @param weight_var Name of the weight variable
#' @param lambda A numeric value for the lambda parameter
#' @param alpha A numeric value for the alpha parameter
#' @param family A character value for the glmnet family argument
#' @param ... Additional arguments to pass to glmnet
#' @return List of objects: glmnet object, data table of betas, ggplot2 object
#' @importFrom dplyr select pull filter
#' @importFrom glmnet glmnet
#' @importFrom data.table copy
#' @export
tidy_glmnet <- function(
    data,
    exposures,
    weight_var = NULL,
    outcome,
    lambda,
    alpha,
    family = "binomial",
    ...
) {
    dataset <- data.table::copy(data)
    dataset <- dataset |>
        (\(x) {
            if (!is.null(weight_var)) {
                x |> dplyr::filter(!is.na(get(weight_var)))
            } else {
                x
            }
        })()
    if (!is.null(weight_var)) {
        weight <- dataset |>
            dplyr::pull(get(weight_var))
    } else {
        weight <- NULL
    }
    model_fit <- glmnet::glmnet(
        x = dataset |>
            dplyr::select(tidyselect::any_of(exposures)) |>
            as.matrix(),
        y            = dataset |> dplyr::pull(outcome),
        weights      = weight,
        alpha        = alpha,
        family       = family,
        lambda       = lambda,
        ...
    )

    model_betas <- tidy_glmnet_betas(model_fit, ...)

    vip_plot <- plot_glmnet_vip(model_betas, exposure_var = label, beta_var = rel_beta)

    return(list(
        model_fit   = model_fit,
        model_betas = model_betas,
        vip_plot    = vip_plot
    ))
}

#' Extract variable importance from ranger model
#' @param x A ranger model object
#' @param phecode_description A logical indicating whether to include the phecode description
#' @return A data table with the exposure and importance
#' @importFrom data.table data.table
#' @importFrom dplyr mutate arrange desc left_join n select
#' @importFrom stringr str_wrap
#' @importFrom ms pheinfox
#' @export
tidy_ranger_imp <- function(x, phecode_description = TRUE) {
        data.table::data.table(
            exposure       = names(x$variable.importance),
            importance     = x$variable.importance
        ) |>    
        dplyr::mutate(
            abs_importance = abs(importance),
            rel_importance = abs_importance / max(abs_importance, na.rm = TRUE),
            sign = ifelse(importance < 0, "Negative", "Positive"),
        ) |>
        dplyr::arrange(dplyr::desc(abs_importance)) |>
        dplyr::mutate(rank = 1:dplyr::n()) |>
        (\(x) if (phecode_description) {
            x |>
                dplyr::left_join(
                    ms::pheinfox |> dplyr::select(exposure = phecode, description),
                    by = "exposure"
                ) |>
                dplyr::mutate(
                    label = stringr::str_wrap(paste0(description, " [", exposure, "]"), width = 35)
                )
        })()
}

#' Plot variable importance from ranger model
#' @param x Table of variable importance from ranger model (see tidy_ranger_imp)
#' @param exposure_var The name of the exposure variable
#' @param top_n The number of top variables to plot
#' @param caption Plot caption
#' @return A ggplot2 object
#' @importFrom dplyr slice_head
#' @importFrom ggplot2 ggplot aes geom_bar labs coord_flip
#' @importFrom stats reorder
#' @export
plot_ranger_vip <- function(x, exposure_var = label, top_n = 15, caption = ggplot2::waiver()) {
    x |>
        dplyr::slice_head(n = top_n) |>
        ggplot2::ggplot(ggplot2::aes(x = stats::reorder(label, rel_importance), y = rel_importance, fill = sign)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
            x       = "Variables",
            y       = "Relative importance",
            title   = paste0("Variable Importance (top ", top_n, ")"),
            caption = caption
        ) +
        ggplot2::theme_minimal() +
        ggplot2::coord_flip()
}

#' Fit a ranger model and extract variable importance and plot
#' @param formula A formula
#' @param data A data frame
#' @param num_trees A numeric value for the number of trees
#' @param mtry A numeric value for the mtry parameter
#' @param min.node.size A numeric value for the min.node.size parameter
#' @param num_threads A numeric value for the num_threads parameter
#' @param importance A character value for the importance parameter
#' @param write_forest A logical value for the write_forest parameter
#' @param exposure_var The name of the exposure variable
#' @param ... Additional arguments to pass to ranger
#' @return List of objects: ranger object, data table of importance, ggplot2 object
#' @importFrom dplyr select pull
#' @importFrom ranger ranger
#' @export
tidy_ranger <- function(
    formula,
    data,
    num_trees = 500,
    mtry,
    min.node.size,
    num_threads  = 12,
    importance   = "permutation",
    write_forest = TRUE,
    exposure_var = label,
    ...
) {
    model_fit <- ranger::ranger(
        formula       = formula,
        data          = data,
        num.trees     = num_trees,
        mtry          = mtry,
        min.node.size = min.node.size,
        num.threads   = num_threads,
        importance    = importance,
        write.forest  = write_forest,
        ...
    )

    model_imp <- tidy_ranger_imp(model_fit)

    vip_plot <- plot_ranger_vip(model_imp, exposure_var = exposure_var, caption = paste0("mtry = ", mtry, ", min.node.size = ", min.node.size, ", importance = ", importance))

    return(list(
        model_fit = model_fit,
        model_imp = model_imp,
        vip_plot  = vip_plot
    ))

}



#' Fit a weighted lasso model and extract variable importance (betas) and plot
#' @param data A data frame
#' @param exposures A character vector of exposure variables
#' @param outcome A character vector of outcome variables
#' @param weight_var Name of the weight variable
#' @param lambda A numeric value for the lambda parameter
#' @param family A character value for the glmnet family argument
#' @param ... Additional arguments to pass to glmnet
#' @return List of objects: glmnet object, data table of betas, ggplot2 object
#' @importFrom dplyr select pull filter
#' @importFrom glmnet glmnet
#' @importFrom data.table copy
#' @importFrom cli cli_alert_danger
#' @export
tidy_wlasso <- function(
    data,
    exposures,
    outcome,
    weight_var = NULL,
    lambda     = NULL,
    family     = "binomial",
    ...
) {
    cols <- c(exposures, outcome, weight_var)
    dataset <- data.table::copy(data)
    dataset <- dataset |>
        (\(x) {
            if (!is.null(weight_var)) {
                x |> dplyr::filter(!is.na(get(weight_var)))
            } else {
                x
            }
        })()

    model_fit <- tryCatch({
        weighted_lasso(
            data        = dataset,
            col.x       = exposures,
            col.y       = outcome,
            weights     = weight_var,
            method      = "dCV",
            lambda.grid = lambda,
            family      = family
        )[["model.min"]]
    },
    error = function(e) {
        cli::cli_alert_danger("Likely convergence issue in wlasso - falling back to glmnet. Error: ", e$message)
        weight <- dataset |>
            dplyr::pull(get(weight_var))
        glmnet::glmnet(
            x = dataset |>
                dplyr::select(tidyselect::any_of(exposures)) |>
                as.matrix(),
            y = dataset |> dplyr::pull(outcome),
            weights = weight,
            alpha = 1,
            family = family,
            lambda = lambda,
            ...
        )
    },
    warning = function(w) {
        cli::cli_alert_danger("Likely convergence issue in wlasso - falling back to glmnet. Warning: ", w$message)
        weight <- dataset |>
            dplyr::pull(get(weight_var))
        glmnet::glmnet(
            x = dataset |>
                dplyr::select(tidyselect::any_of(exposures)) |>
                as.matrix(),
            y = dataset |> dplyr::pull(outcome),
            weights = weight,
            alpha = 1,
            family = family,
            lambda = lambda,
            ...
        )
    })

    model_betas <- tidy_glmnet_betas(model_fit, ...)

    vip_plot <- plot_glmnet_vip(model_betas, exposure_var = label, beta_var = rel_beta)

    return(list(
        model_fit   = model_fit,
        model_betas = model_betas,
        vip_plot    = vip_plot
    ))
}
