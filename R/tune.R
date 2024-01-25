#' Tune lambda (and alpha) hyperparameters for lasso, ridge, and elastic net regression via glmnet
#' @param data dataset
#' @param outcome name of outcome variable
#' @param exposures vector of exposure variable names; if NULL, assumes all non-outcome variables
#' @param alpha 0 for lasso, 1 for ridge, or a vector of numeric values to try for elastic net
#' @param parallel parallelize?
#' @param n_folds number of folds for cross-validation
#' @param family model family (binomial or guassian)
#' @param verbose print additional information
#' @param return_mod return the full cv.glmnet object?
#' @param ... additional arguments passed to cv.glmnet
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table
#' @importFrom cli cli_progress_step cli_progress_done
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @return return a table with hyperparameters and their values
#' @export
tune_glmnet <- function(
    data,
    outcome,
    exposures,
    .alpha,
    parallel = TRUE,
    n_folds = 10,
    family  = "binomial",
    verbose = TRUE,
    return_mod = FALSE,
    ...
) {

  if (!data.table::is.data.table(data)) data <- data.table::as.data.table(data)
  x <- as.matrix(data |> dplyr::select(tidyselect::all_of(exposures)))
  y <- data[[outcome]]

  set.seed(123)

  if (length(.alpha) == 1) {
    if (.alpha == 1) {
      if (verbose) cli::cli_progress_step("Fitting lasso...")
      prefix <- "lasso_"
    }
    if (.alpha == 0) {
      if (verbose) cli::cli_progress_step("Fitting ridge...")
      prefix <- "ridge_"
    }
    .alpha_grid = .alpha
  } else {
    if (verbose) cli::cli_progress_step("Fitting elastic net...")
    .alpha_grid <- .alpha
    prefix     <- "enet_"
  }
  if (verbose) on.exit(cli::cli_progress_done())

  cv_fit_list <- list()

  for (a in .alpha_grid) {
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

  out <- list(
    param = data.table::data.table(
      "parameter" = paste0(prefix, c("alpha", "lambda.min", "lambda.1se")),
      "value"     = c(best_alpha, best_lambda, best_lambda_1se)
    )
  )

  if (return_mod) {
    out$mod <- glmnet::glmnet(
      x, y,
      alpha = best_alpha,
      lambda = best_lambda,
      nfolds = n_folds,
      family = family,
      parallel = parallel,
    )
  }

  return(out)

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
#' @param return_mod return the full ranger object?
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
    verbose       = TRUE,
    return_mod   = FALSE
) {
  if (verbose) {
    cli::cli_progress_step("Fitting ranger random forest...")
    on.exit(cli::cli_progress_done())
  }

  if (is.null(mtry_seq)) {
    p <- ncol(data) - 1
    mtry_seq <- floor(sqrt(p) * c(0.25, 0.5, 1, 2, 4))
  }

  if (is.null(node_size_seq)) {
    min.node.size_grid <- c(1, 3, 5, 10, 20)
    min.node.size_grid <- min.node.size_grid[min.node.size_grid < nrow(data)]
  }

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

  out <- list(
    param = data.table::data.table(
      "parameter" = c("rf.mtry", "rf.node_size"),
      "value" = c(
        rf |>
          dplyr::filter(oob_rmse == min(oob_rmse, na.rm = TRUE)) |>
          dplyr::pull(mtry),
        rf |>
          dplyr::filter(oob_rmse == min(oob_rmse, na.rm = TRUE)) |>
          dplyr::pull(node_size)
      )
    )
  )

  if (return_mod) {
    model <- ranger::ranger(
      formula         = f,
      data            = data |> dplyr::select(tidyselect::all_of(vars)),
      num.trees       = n_trees,
      mtry            = out$param[out$param$parameter == "rf.mtry", ][["value"]],
      min.node.size   = out$param[out$param$parameter == "rf.node_size", ][["value"]],
      num.threads     = n_cores,
      write.forest = TRUE,
      importance = "permutation"
    )
    out$mod <- model
  }

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
#' @param return_mod return the full wglmnet object?
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
    parallel = FALSE,
    verbose = TRUE,
    method = "dCV",
    family = "binomial",
    return_mod = FALSE,
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
    out <- list(
      param = data.table::data.table(
        "parameter" = paste0(prefix, c("lambda.min", "alpha")),
        "value"     = c(wglmnet_fit$lambda$min, wglmnet_fit$alpha$min)
      )
    )
    if (return_mod) {
      out$mod <- wglmnet_fit
    }
    return(out)
  }, error = function(err_msg) {
    print(err_msg)
    cli::cli_alert_danger("issue with wglmnet, switching to glmnet")
    ex_check  <- exposures[!(exposures %in% outcome)]
    wex_check <- c(ex_check, weight)
    pf        <- as.numeric(wex_check != weight)
    tune_glmnet(
      data           = data,
      outcome        = outcome,
      exposures      = wex_check,
      n_folds        = n_folds,
      .alpha         = a,
      penalty.factor = pf,
      parallel       = parallel,
      return_mod    = return_mod,
    )
  }, warning = function(wrn_msg) {
    cli::cli_alert_danger("issue with wglmnet, switching to glmnet: {wrn_msg}")
    ex_check  <- exposures[!(exposures %in% outcome)]
    wex_check <- c(ex_check, weight)
    pf        <- as.numeric(wex_check != weight)
    tune_glmnet(
      data           = data,
      outcome        = outcome,
      exposures      = wex_check,
      n_folds        = n_folds,
      .alpha         = a,
      penalty.factor = pf,
      return_mod = return_mod,
      parallel       = parallel
    )
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
#' @param return_mod return the full model object?
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
    verbose       = TRUE,
    return_mod    = TRUE
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
    use_these_vars <- names(data)[!names(data) %in% c("id", weight)]
    dataset  <- data.table::copy(data[stats::complete.cases(data[use_these_vars]), ])
    dataset[["id"]] <- NULL
    wdataset <- data.table::copy(dataset[stats::complete.cases(dataset), ])
    if (is.null(exposures)) {
      exposures <- names(dataset)[!names(dataset) %in% c(outcome, weight)]
    }
  }

  out <- list()

  if ("ridge" %in% methods) {
    ridge_mod <- tune_glmnet(
      data      = dataset,
      outcome   = outcome,
      exposures = exposures,
      n_folds   = n_folds,
      .alpha    = 0,
      parallel  = parallel,
      verbose   = verbose,
      return_mod = return_mod
    )
    out[["ridge"]] <- ridge_mod
    if (!is.null(weight)) {
      wridge_mod <- tune_wglmnet(
        data           = wdataset,
        outcome        = outcome,
        exposures      = exposures,
        weight         = weight,
        n_folds        = n_folds,
        alpha          = 0,
        parallel       = parallel,
        verbose        = verbose,
      return_mod = return_mod
      ) 
      out[["wridge"]] <- wridge_mod
    }
  }

  if ("lasso" %in% methods) {
    lasso_mod <- tune_glmnet(
      data      = dataset,
      outcome   = outcome,
      exposures = exposures,
      n_folds   = n_folds,
      .alpha    = 1,
      parallel  = parallel,
      verbose   = verbose,
      return_mod = return_mod
    )
    out[["lasso"]] <- lasso_mod
    if (!is.null(weight)) {
      wlasso_mod <- tune_wglmnet(
        data      = wdataset,
        outcome   = outcome,
        exposures = exposures,
        weight    = weight,
        alpha     = 1,
        n_folds   = n_folds,
        verbose   = verbose,
      return_mod = return_mod
      )
      out[["wlasso"]] <- wlasso_mod
    }
  }

  if ("enet" %in% methods) {
    if (is.null(alpha)) alpha <- seq(0, 1, length.out = 20)
    enet_mod <- tune_glmnet(
      data      = dataset,
      outcome   = outcome,
      exposures = exposures,
      n_folds   = n_folds,
      .alpha    = alpha,
      parallel  = parallel,
      verbose   = verbose,
      return_mod = return_mod
    )
    out[["enet"]] <- enet_mod
    if (!is.null(weight)) {
      wenet_mod <- tune_wglmnet(
        data           = wdataset,
        outcome        = outcome,
        exposures      = exposures,
        weight         = weight,
        n_folds        = n_folds,
        alpha          = alpha,
        parallel       = parallel,
        verbose        = verbose,
      return_mod = return_mod
      )
      out[["wenet"]] <- wenet_mod
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
      verbose       = verbose,
      return_mod = return_mod
    )
    out[["rf"]] <- rf_mod
    if (!is.null(weight)) {
      wrf_mod <- wglmnet::wranger(
        data = wdataset,
        col.y = outcome,
        col.x = exposures,
        weights = weight,
        mtry_grid             = mtry_seq,
    min.node.size_grid = node_size_seq,
    num_trees             = n_trees
      )
      out[["wrf"]] <- wrf_mod
    }
  }

  return(out)

}
