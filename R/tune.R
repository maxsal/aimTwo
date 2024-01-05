#' Tune lambda (and alpha) hyperparameters for lasso, ridge, and elastic net regression via glmnet
#' @param data dataset
#' @param outcome name of outcome variable
#' @param exposures vector of exposure variable names; if NULL, assumes all non-outcome variables
#' @param alpha 0 for lasso, 1 for ridge, or a vector of numeric values to try for elastic net
#' @param parallel parallelize?
#' @param n_folds number of folds for cross-validation
#' @param family model family (binomial or guassian)
#' @param verbose print additional information
#' @param ... additional arguments passed to cv.glmnet
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table
#' @importFrom cli cli_progress_step cli_progress_done
#' @return return a table with hyperparameters and their values
#' @export
tune_glmnet <- function(
    data,
    outcome,
    exposures,
    alpha,
    parallel = TRUE,
    n_folds = 10,
    family  = "binomial",
    verbose = TRUE,
    ...
) {

  x <- as.matrix(data[exposures])
  y <- data[[outcome]]

  set.seed(123)

  if (length(alpha) == 1) {
    if (alpha == 1) {
      if (verbose) cli::cli_progress_step("Fitting lasso...")
      prefix <- "lasso_"
    }
    if (alpha == 0) {
      if (verbose) cli::cli_progress_step("Fitting ridge...")
      prefix <- "ridge_"
    }
    alpha_grid = alpha
  } else {
    if (verbose) cli::cli_progress_step("Fitting elastic net...")
    alpha_grid <- alpha
    prefix     <- "enet_"
  }
  if (verbose) on.exit(cli::cli_progress_done())

  cv_fit_list <- list()

  for (a in alpha_grid) {
    cv_fit <- glmnet::cv.glmnet(
      x, y,
      alpha    = a,
      nfolds   = n_folds,
      family   = family,
      parallel = parallel,
      ...
    )

    cv_fit$alpha <- a
    cv_fit_list[[paste0("alpha_", a)]] <- cv_fit
  }

  # Find model with smallest cvm
  min_cvm         <- min(sapply(cv_fit_list, function(x) min(x$cvm)))
  best_model      <- cv_fit_list[[which.min(sapply(cv_fit_list, function(x) min(x$cvm)))]]
  best_alpha      <- best_model$alpha
  best_lambda     <- best_model$lambda.min
  best_lambda_1se <- best_model$lambda.1se

  return(data.table::data.table(
    "parameter" = paste0(prefix, c("alpha", "lambda.min", "lambda.1se")),
    "value"     = c(best_alpha, best_lambda, best_lambda_1se)
  ))
}


#' Tune mtry and min.node.size parameters for random forest via ranger
#' @param data dataset
#' @param outcome name of outcome variable
#' @param exposures vector of exposure variable names; if NULL, assumes all non-outcome variables
#' @param mtry_seq values of mtry to search over. defaults to some recommendations from Breiman
#' @param node_size_seq values of min.node.size to search over. defaults to 5 values between 3 and nrow(data)/100
#' @param n_cores number of cores if parallel
#' @param n_trees number of trees to grow
#' @param verbose print additional information
#' @importFrom ranger ranger
#' @importFrom data.table data.table rbindlist as.data.table
#' @importFrom cli  cli_progress_step cli_alert cli_progress_update cli_progress_done
#' @importFrom dplyr select filter pull
#' @importFrom tidyselect all_of
#' @return return a table with hyperparameters and their values
#' @export
tune_ranger <- function(
    data,
    outcome,
    exposures     = NULL,
    mtry_seq      = NULL,
    node_size_seq = NULL,
    n_cores       = 1,
    n_trees       = 500,
    verbose       = TRUE
) {
  if (verbose) {
    cli::cli_progress_step("Fitting ranger random forest...")
    on.exit(cli::cli_progress_done())
  }
  out <- data.table::data.table()

  if (is.null(mtry_seq)) {
    p <- ncol(data) - 1
    mtry_seq <- round(c(p/2, p/3, sqrt(p) * c(0.5, 1, 2), log2(p)))
  }

  if (is.null(node_size_seq)) {
    node_size_seq <- round(seq(3, round(nrow(data) / 100), length.out = 5))
  }

  cli::cli_progress_step("Fitting random forest...")
  rf <- expand.grid(
    mtry        = mtry_seq,
    node_size   = node_size_seq,
    oob_rmse    = NA
  ) |> data.table::as.data.table()

  cli::cli_alert("Trying {length(mtry_seq)} values of mtry, {length(node_size_seq)} of node size ({nrow(rf)} combinations)")

  if (is.null(exposures)) exposures <- names(data)[names(data) != outcome]
  f <- formula(paste0(outcome, " ~ ", paste0(exposures, collapse = " + ")))
  vars <- c(outcome, exposures)

  cli::cli_progress_bar(name = "grid search...", total = nrow(rf))
  for (i in 1:nrow(rf)) {
    model <- ranger::ranger(
      formula         = f,
      data            = data |> dplyr::select(tidyselect::all_of(vars)),
      num.trees       = n_trees,
      mtry            = rf$mtry[i],
      min.node.size   = rf$node_size[i],
      num.threads     = n_cores
    )
    rf$oob_rmse[i] <- sqrt(model$prediction.error)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  out <- data.table::rbindlist(list(
    out,
    data.table::data.table(
      "parameter" = c("rf.mtry", "rf.node_size"),
      "value"     = c(
        rf |>
            dplyr::filter(oob_rmse == min(oob_rmse, na.rm = TRUE)) |>
            dplyr::pull(mtry),
        rf |>
            dplyr::filter(oob_rmse == min(oob_rmse, na.rm = TRUE)) |>
            dplyr::pull(node_size)
      )
    )
  ), use.names = TRUE, fill = TRUE)

  return(out)

}

#' Tune lambda parameter for weighted lasso via wlasso::wlasso (then fall back to glmnet if issues)
#' @param data dataset
#' @param outcome name of outcome variable
#' @param exposures vector of exposure variable names; if NULL, assumes all non-outcome variables
#' @param alpha 0 for lasso, 1 for ridge, or a vector of numeric values to try for elastic net
#' @param weight name of weight variable
#' @param n_folds number of folds for cross-validation
#' @param parallel parallelize?
#' @param verbose print additional information
#' @param method method for replicate weight sampling
#' @param family model family (binomial or guassian)
#' @param ... additional arguments passed to wglmnet
#' @importFrom wglmnet wglmnet
#' @importFrom cli cli_alert_danger cli_progress_step cli_progress_done
#' @importFrom dplyr mutate
#' @importFrom data.table data.table
#' @return return a table with hyperparameters and their values
#' @export
tune_wglmnet <- function(
    data,
    outcome,
    exposures,
    alpha = 1,
    weight,
    n_folds = 10,
    parallel = parallel,
    verbose = TRUE,
    method = "dCV",
    family = "binomial",
    ...
) {
  if (length(alpha) == 1) {
    if (alpha == 1) {
      if (verbose) cli::cli_progress_step("Fitting weighted lasso...")
      prefix <- "wlasso_"
    }
    if (alpha == 0) {
      if (verbose) cli::cli_progress_step("Fitting weighted ridge...")
      prefix <- "wridge_"
    }
    alpha_grid <- alpha
  } else {
    if (verbose) cli::cli_progress_step("Fitting weighted elastic net...")
    alpha_grid <- alpha
    prefix <- "wenet_"
  }
  if (verbose) on.exit(cli::cli_progress_done())
  a <- alpha

  tryCatch({
    wglmnet_fit <- wglmnet::wglmnet(
      data = data, col.y = outcome,
      col.x = exposures[!(exposures %in% c(outcome, weight))],
      alpha = alpha,
      weights = weight,
      family = family,
      lambda.grid = NULL,
      method = method,
      k = n_folds,
      ...
    )
    data.table::data.table(
      "parameter" = paste0(prefix, c("lambda.min", "alpha")),
      "value"     = c(wglmnet_fit$lambda$min, wglmnet_fit$alpha$min)
    )
  }, error = function(err_msg) {
    print(err_msg)
    cli::cli_alert_danger("issue with wglmnet, switching to glmnet")
    pf <- as.numeric(names(data)[names(data) != outcome] != weight)
    tune_glmnet(
      data           = data,
      outcome        = outcome,
      exposures      = exposures[!(exposures %in% outcome)],
      n_folds        = n_folds,
      alpha          = a,
      penalty_factor = pf,
      parallel       = parallel
    ) |> dplyr::mutate(parameter = paste0("w", parameter))
  }, warning = function(wrn_msg) {
    cli::cli_alert_danger("issue with wglmnet, switching to glmnet: {wrn_msg}")
    pf <- as.numeric(names(data)[names(data) != outcome] != weight)
    tune_glmnet(
      data           = data,
      outcome        = outcome,
      exposures      = exposures[!(exposures %in% outcome)],
      n_folds        = n_folds,
      alpha          = a,
      penalty_factor = pf,
      parallel       = parallel
    ) |> dplyr::mutate(parameter = paste0("w", parameter))
  })
}

#' Tune hyperparameters for (un)weighted regularized regression and unweighted random forest
#' @param data dataset
#' @param outcome name of outcome variable
#' @param exposures vector of exposure variable names; if NULL, assumes all non-outcome variables
#' @param weight name of weight variable
#' @param parallel parallelize?
#' @param methods vector of method names to consider: ridge, lasso, enet, rf
#' @param n_folds number of folds for cross-validation
#' @param n_cores number of cores to use for parallelization
#' @param alpha value(s) of alpha to pass to tune_glmnet if ridge, lasso, and/or enet are selected
#' @param mtry_seq values of mtry to try via tune_rf if rf method is selected
#' @param node_size_seq value of min.node.size to try via tune_rf if rf method is selected
#' @param n_trees number of trees to grow if rf method is selected
#' @param verbose print additional information
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr mutate bind_rows
#' @importFrom data.table copy data.table
#' @importFrom stats complete.cases
#' @importFrom cli cli_progress_step cli_progress_done
#' @return return a table with hyperparameters and their values
#' @export
tune_models <- function(
    data,
    outcome,
    exposures     = NULL,
    weight        = NULL,
    parallel      = TRUE,
    methods       = c("ridge", "lasso", "enet", "rf"),
    n_folds       = 10,
    n_cores       = 1,
    alpha         = NULL,
    mtry_seq      = NULL,
    node_size_seq = NULL,
    n_trees       = 500,
    verbose       = TRUE
) {

  if (parallel) {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }

  if (is.null(weight)) {
    dataset <- data.table::copy(data[stats::complete.cases(data), ])
    if (is.null(exposures)) {
      exposures <- names(dataset)[names(dataset) != outcome]
    }
  } else {
    dataset  <- data.table::copy(data[stats::complete.cases(data[names(data)[names(data) != weight]]), ])
    wdataset <- data.table::copy(dataset[stats::complete.cases(dataset), ])
    if (is.null(exposures)) {
      exposures <- names(dataset)[!names(dataset) %in% c(outcome, weight)]
    }
  }

  out <- data.table::data.table()

  if ("ridge" %in% methods) {
    ridge_mod <- tune_glmnet(
      data      = dataset,
      outcome   = outcome,
      exposures = exposures,
      n_folds   = n_folds,
      alpha     = 0,
      parallel  = parallel,
      verbose   = verbose
    )
    out <- dplyr::bind_rows(out, ridge_mod)
    if (!is.null(weight)) {
      wridge_mod <- tune_wglmnet(
        data           = wdataset,
        outcome        = outcome,
        exposures      = exposures,
        weight         = weight,
        n_folds        = n_folds,
        alpha          = 0,
        parallel       = parallel,
        verbose        = verbose
      ) |>
        dplyr::mutate(parameter = paste0("w", parameter))
      out <- dplyr::bind_rows(out, wridge_mod)
    }
  }

  if ("lasso" %in% methods) {
    lasso_mod <- tune_glmnet(
      data      = dataset,
      outcome   = outcome,
      exposures = exposures,
      n_folds   = n_folds,
      alpha     = 1,
      parallel  = parallel,
      verbose   = verbose
    )
    out <- dplyr::bind_rows(out, lasso_mod)
    if (!is.null(weight)) {
      wlasso_mod <- tune_wglmnet(
        data      = wdataset,
        outcome   = outcome,
        exposures = exposures,
        weight    = weight,
        alpha     = 1,
        n_folds   = n_folds,
        verbose   = verbose
      )
      out <- dplyr::bind_rows(out, wlasso_mod)
    }
  }

  if ("enet" %in% methods) {
    if (is.null(alpha)) alpha <- seq(0, 1, length.out = 20)
    enet_mod <- tune_glmnet(
      data      = dataset,
      outcome   = outcome,
      exposures = exposures,
      n_folds   = n_folds,
      alpha     = alpha,
      parallel  = parallel,
      verbose   = verbose
    )
    out <- dplyr::bind_rows(out, enet_mod)
    if (!is.null(weight)) {
      wenet_mod <- tune_wglmnet(
        data           = wdataset,
        outcome        = outcome,
        exposures      = exposures,
        weight         = weight,
        n_folds        = n_folds,
        alpha          = alpha,
        parallel       = parallel,
        verbose        = verbose
      ) |> dplyr::mutate(parameter = paste0("w", parameter))
      out <- dplyr::bind_rows(out, wenet_mod)
    }
  }

  if ("rf" %in% methods) {
    rf_mod <- tune_ranger(
      data          = dataset,
      outcome       = outcome,
      exposures     = exposures,
      n_cores       = n_cores,
      mtry_seq      = mtry_seq,
      node_size_seq = node_size_seq,
      n_trees       = n_trees,
      verbose       = verbose
    )
    out <- dplyr::bind_rows(out, rf_mod)
  }

  return(out)

}
