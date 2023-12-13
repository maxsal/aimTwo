#' Replicate weights (based on wlasso::replicate.weights())
#'
#' @description This function allows calculating replicate weights.
#'
#' @param data A data frame with information of cluster and strata indicators and sampling weights.
#' @param method To choose between \code{JKn}, \code{dCV}, \code{bootstrap}, \code{subbootstrap}, \code{BRR}, \code{split} or \code{extrapolation}
#' @param cluster Name of the column indicating clusters
#' @param strata Name of the column indicating strata
#' @param weights Name of the column indicating sampling weights
#' @param k Number of folds, if \code{dCV} method is selected
#' @param B Number of bootstrap resamples if \code{bootstrap} or \code{subbootstrap} methods are selected
#' @param R Number of times the sample is split (only for \code{cv}, \code{split} and \code{extrapolation} methods)
#' @param train_prob Proportion of PSUs set as training set.
#' @param method_split To choose between \code{dCV}, \code{bootstrap} or \code{subbootstrap}
#' @param seed seed
#'
#' @importFrom stats as.formula
#' @importFrom survey svydesign as.svrepdesign
#' 
#' @return This function returns a new data frame with new columns, each of them indicating replicate weights for different subsets.
#' @export

replicate_weights <- function(
  data,
  method       = "bootstrap",
  cluster      = NULL,
  strata       = NULL,
  weights      = NULL,
  k            = 10,
  B            = 200,
  R            = 1,
  train_prob   = NULL,
  method_split = c("dCV", "bootstrap", "subbootstrap"),
  seed         = NULL
) {

  if (method == "dCV") {

    new_data <- cv_folds(data, k, weights, seed, strata, cluster, R)

    for(r in 1:R){
      for(kk in 1:k) {
        new_data[, paste0("sw_r_", r, "_test_", kk)] <- rep(0, nrow(new_data))
        new_data[which(new_data[,paste0("folds_", r)] == kk),
                       paste0("sw_r_", r, "_test_", kk)] <- new_data[which(new_data[, paste0("folds_", r)] == kk), weights]

      }
    }

  } else {

    if (method == "split") {
      new_data <- rw_split(data, train_prob, method = method_split,
                          weights, strata, cluster, R, seed)
    } else {

      if (method == "extrapolation") {
        new_data <- split_strata(data, train_prob, strata, weights, seed, R)
      } else {

      # Define cluster formula
      if (is.null(cluster)) {
        formula_cluster <- stats::as.formula("~1")
      } else {
        formula_cluster <- stats::as.formula(paste0("~", cluster))
      }

      # Define strata formula
      if (!is.null(strata)) {
        formula_strata <- stats::as.formula(paste0("~", strata))
      } else {
        formula_strata <- NULL
      }

      # Define weights formula
      if (!is.null(weights)) {
        formula_weights <- stats::as.formula(paste0("~", weights))
      }

      # Define the design
      des <- survey::svydesign(
        ids     = formula_cluster,
        strata  = formula_strata,
        weights = formula_weights,
        data    = data,
        nest    = FALSE
      )

      # Generate replicate weights based on the selected method
      if (method %in% c("JKn", "bootstrap", "subbootstrap", "BRR")) {
        if (method %in% c("bootstrap", "subbootstrap")) {
          rep_des <- survey::as.svrepdesign(design = des, type = method, replicates = B)
        } else {
          if (is.null(strata) & method == "JKn") {
            rep_des <- survey::as.svrepdesign(design = des, type = "JK1")
          } else {
            rep_des <- survey::as.svrepdesign(design = des, type = method)
          }
        }

        mat_repw_ind       <- apply(rep_des$repweights$weights, 2, function(x) { x[rep_des$repweights$index] })
        mat_repw           <- apply(mat_repw_ind, 2, function(x) { x * data[, weights] } )
        colnames(mat_repw) <- paste0("rw_r_1_train_", 1:ncol(mat_repw))
        new_data           <- cbind(data, mat_repw)
      }

      # Define replicate weights for the testing set in BRR
      if (method == "BRR") {
        mat_repw_test           <- mat_repw_ind
        colnames(mat_repw_test) <- paste0("rw_r_1_test_", 1:ncol(mat_repw_test))
        mat_repw_test           <- -1 * (mat_repw_test - 2)
        mat_repw_test           <- apply(mat_repw_test, 2, function(x) { x * data[, weights] })
        new_data                <- cbind(new_data, mat_repw_test)
      }

      # Define each unit as fold in JKn method
      if (method == "JKn") {
        mat_repw_test                           <- -mat_repw_ind
        mat_repw_test[which(mat_repw_test==0)]  <- 1
        mat_repw_test[which(mat_repw_test < 0)] <- 0
        mat_repw_test                           <- apply(mat_repw_test, 2, function(x){ x * data[,weights] })
        colnames(mat_repw_test)                 <- paste0("sw_r_1_test_", 1:ncol(mat_repw_test))
        new_data                                <- cbind(new_data, mat_repw_test)
      }
    }
    }
  }
  return(new_data)
}



# dCV
#' Cross-validation folds
#' @description This function allows calculating cross-validation folds.
#' @param data A data frame with information of cluster and strata indicators and sampling weights.
#' @param k Number of folds
#' @param weights Name of the column indicating sampling weights
#' @param seed seed
#' @param strata Name of the column indicating strata
#' @param cluster Name of the column indicating clusters
#' @param R Number of times the sample is split
#' @return This function returns a new data frame with new columns, each of them indicating the fold of each unit.
#' @export
cv_folds <- function(
  data,
  k,
  weights,
  seed    = NULL,
  strata  = NULL,
  cluster = NULL,
  R       = 1
) {

  set.seed(seed)
  seeds <- stats::runif(R) * 10000

  if (is.null(cluster) & !is.null(strata)) {
    data$cluster <- 1:nrow(data)
    cluster <- "cluster"
  } else {
    if (!is.null(cluster) & is.null(strata)) {
      data$strata <- rep(1, nrow(data))
      strata      <- "strata"
    } else {
      if (is.null(cluster) & is.null(strata)) {
        data$strata  <- rep(1, nrow(data))
        data$cluster <- 1:nrow(data)
        strata       <- "strata"
        cluster      <- "cluster"
      }
    }
  }

  for (r in 1:R) {
    data[,paste0("folds_",r)] <- f_folds(
      data,
      k       = k,
      seed    = seeds[r],
      strata  = strata,
      cluster = cluster
    )
    for (kk in 1:k) {
      data[, paste0("rw_r_", r, "_train_", kk)] <- repl_weights(
        data,
        folds     = paste0("folds_", r),
        test_fold = kk,
        weights,
        strata,
        cluster)
    }

    for (kk in 1:k) {
      data[,paste0("rw_r_",r,"_test_", kk)] <- repl_weights_test(
        data,
        folds     = paste0("folds_", r),
        test_fold = kk,
        weights,
        strata,
        cluster
      )
    }
  }
  return(data)
}

#' Folds for cross-validation
#' @description This function allows calculating folds for cross-validation.
#' @param data A data frame with information of cluster and strata indicators and sampling weights.
#' @param k Number of folds
#' @param seed seed
#' @param strata Name of the column indicating strata
#' @param cluster Name of the column indicating clusters
#' @return This function returns a vector with the fold of each unit.
#' @export
f_folds <- function(
  data,
  k       = 5,
  strata  = NULL,
  cluster = NULL,
  seed    = NULL
) {
  if(!is.null(seed)){ set.seed(seed) }

  data$hclus <- interaction(data[,strata], data[,cluster], drop = TRUE)
  new_ids    <- levels(data$hclus)
  n          <- length(new_ids)

  v_folds_new_ids        <- sample(cut(seq(1, n), breaks = k, labels = FALSE))
  names(v_folds_new_ids) <- new_ids
  v_folds                <- v_folds_new_ids[match(data$hclus, names(v_folds_new_ids))]

  h_fold     <- table(data[, strata], v_folds) != 0
  h_fold_sum <- apply(h_fold, 1, sum)
  h_onefold  <- which(h_fold_sum == 1)
  if (length(h_onefold) != 0) {
    for (hh in h_onefold) {
      kk                       <- which(h_fold[hh,] == 1)
      id_h                     <- which(data[, strata] == hh)
      psu_h                    <- unique(data$hclus[id_h])
      set.seed(seed * 10)
      selected_psu             <- sample(psu_h, 1)
      set.seed(seed * 100)
      new_k                    <- sample(c(1:k)[-kk], 1)
      id_selected_psu          <- which(data$hclus == selected_psu)
      v_folds[id_selected_psu] <- new_k
    }
  }
  return(v_folds)
}

#' Replicate weights for training set
#' @description This function allows calculating replicate weights for training set.
#' @param data A data frame with information of cluster and strata indicators and sampling weights.
#' @param folds Name of the column indicating the fold of each unit.
#' @param test_fold Fold of the testing set.
#' @param weights Name of the column indicating sampling weights
#' @param strata Name of the column indicating strata
#' @param cluster Name of the column indicating clusters
#' @return This function returns a vector with the replicate weights for the training set.
#' @export
repl_weights <- function(
  data,
  folds,
  test_fold,
  weights,
  strata = NULL,
  cluster = NULL
) {
  v_repl_weights <- rep(0, nrow(data))

  id_test <- which(data[, folds] == test_fold)

  data[, strata] <- as.factor(data[, strata])

  str_clus      <- table(data[, strata], data[, cluster]) != 0
  str_clus_test <- table(data[id_test, strata], data[id_test, cluster]) != 0

  v_mh <- apply(str_clus_test, 1, sum)
  v_nh <- apply(str_clus, 1, sum)
  coef <- v_nh / (v_nh - v_mh)

  v_repl_weights[-id_test] <- data[-id_test, weights] * coef[match(data[-id_test, strata], names(coef))]

  return(v_repl_weights)
}

#' Replicate weights for testing set
#' @description This function allows calculating replicate weights for testing set.
#' @param data A data frame with information of cluster and strata indicators and sampling weights.
#' @param folds Name of the column indicating the fold of each unit.
#' @param test_fold Fold of the testing set.
#' @param weights Name of the column indicating sampling weights
#' @param strata Name of the column indicating strata
#' @param cluster Name of the column indicating clusters
#' @return This function returns a vector with the replicate weights for the testing set.
#' @export
repl_weights_test <- function(
  data,
  folds,
  test_fold,
  weights,
  strata = NULL,
  cluster = NULL
) {
  v_repl_weights <- rep(0, nrow(data))

  id_test <- which(data[,folds]!=test_fold)

  data[, strata] <- as.factor(data[, strata])

  str_clus      <- table(data[, strata], data[, cluster]) != 0
  str_clus_test <- table(data[id_test, strata], data[id_test, cluster]) != 0

  v_mh <- apply(str_clus_test, 1, sum)
  v_nh <- apply(str_clus, 1, sum)
  coef <- v_nh/(v_nh - v_mh)

  v_repl_weights[-id_test] <- data[-id_test,weights]*coef[match(data[-id_test,strata], names(coef))]

  return(v_repl_weights)
}

# Split-sample ------------------------------------------------------------

#' Split-sample
#' @description This function allows calculating split-sample.
#' @param data A data frame with information of cluster and strata indicators and sampling weights.
#' @param train_prob Proportion of PSUs set as training set.
#' @param method To choose between \code{dCV}, \code{bootstrap} or \code{subbootstrap}
#' @param weights Name of the column indicating sampling weights
#' @param strata Name of the column indicating strata
#' @param cluster Name of the column indicating clusters
#' @param R Number of times the sample is split
#' @param seed seed
#' @return This function returns a new data frame with new columns, each of them indicating replicate weights for different subsets.
#' @export
rw_split <- function(
  data,
  train_prob,
  method  = c("dCV", "bootstrap", "subbootstrap"),
  weights,
  strata  = NULL,
  cluster = NULL,
  R       = 1,
  seed    = 1
) {

  set.seed(seed)
  seeds <- stats::runif(R) * 100000

  if (is.null(cluster) & !is.null(strata)) {
    data$cluster <- 1:nrow(data)
    cluster      <- "cluster"
  } else {
    if (!is.null(cluster) & is.null(strata)) {
      data$strata <- rep(1, nrow(data))
      strata      <- "strata"
    } else {
      if (is.null(cluster) & is.null(strata)) {
        data$strata  <- rep(1, nrow(data))
        data$cluster <- 1:nrow(data)
        strata       <- "strata"
        cluster      <- "cluster"
      }
    }
  }

  for(r in 1:R){
    data <- split_sample(data, train_prob, r, strata, cluster, seeds[r])
    tags <- as.vector(unique(data[,paste0("set_r_",r)]))
    if (method == "dCV") {
      for (tag in tags) {
        data[, paste0("rw_r_", r, "_", tag)] <- repl_weights_test(
          data,
          folds     = paste0("set_r_", r),
          test_fold = tag,
          weights,
          strata,
          cluster
        )
      }
    } else {
      if (method %in% c("bootstrap", "subbootstrap")) {
        for (tag in tags) {
          data <- replicate_sample(
            data,
            set = paste0("set_r_", r),
            tag,
            strata,
            weights,
            r,
            boot_type = method
          )
        }
      }
    }
  }

  return(data)
}

#' Split-sample
#' @description This function allows calculating split-sample.
#' @param data A data frame with information of cluster and strata indicators and sampling weights.
#' @param train_prob Proportion of PSUs set as training set.
#' @param r Number of times the sample is split
#' @param strata Name of the column indicating strata
#' @param cluster Name of the column indicating clusters
#' @param seed seed
#' @return This function returns a new data frame with new columns, each of them indicating the fold of each unit.
#' @export
split_sample <- function(
  data,
  train_prob,
  r,
  strata  = NULL,
  cluster = NULL,
  seed    = 1
) {
  data[, strata] <- as.factor(data[, strata])

  set <- paste0("set_r_", r)

  set.seed(seed)

  data$hclus <- interaction(data[, strata], data[, cluster], drop = TRUE)
  new_ids    <- levels(data$hclus)
  n          <- length(new_ids)

  factor             <- c(0, train_prob, 1)
  set_new_ids        <- sample(cut(seq(1,n)/n, factor, labels = c("train", "test")))
  names(set_new_ids) <- new_ids
  data[, set]        <- as.factor(set_new_ids[match(data$hclus, names(set_new_ids))])

  train_0 <- table(data[which(data[, set] == "train"), strata]) == 0
  if (sum(train_0) != 0) {
    h_0 <- which(train_0 == 1)
    for (hh in h_0) {
      id_hh                      <- which(data[, strata] == hh)
      psu_h                      <- unique(data$hclus[id_hh])
      selected_psu               <- sample(psu_h, size = 1)
      id_selected_psu            <- which(data$hclus == selected_psu)
      data[id_selected_psu, set] <- "train"
    }
  }

  test_0 <- table(data[which(data[, set] == "test"), strata]) == 0
  if (sum(test_0) != 0) {
    h_0 <- which(test_0 == 1)
    for (hh in h_0) {
      id_hh                     <- which(data[, strata] == hh)
      psu_h                     <- unique(data$hclus[id_hh])
      selected_psu              <- sample(psu_h, size = 1)
      id_selected_psu           <- which(data$hclus == selected_psu)
      data[id_selected_psu,set] <- "test"
    }
  }

  return(data)
}

#' Replicate sample
#' @description This function allows calculating replicate sample.
#' @param data A data frame with information of cluster and strata indicators and sampling weights.
#' @param set Name of the column indicating the fold of each unit.
#' @param tag Fold of the testing set.
#' @param strata Name of the column indicating strata
#' @param weights Name of the column indicating sampling weights
#' @param r Number of times the sample is split
#' @param boot_type To choose between \code{bootstrap} or \code{subbootstrap}
#' @return This function returns a new data frame with new columns, each of them indicating replicate weights for different subsets.
#' @export
replicate_sample <- function(
  data,
  set,
  tag,
  strata,
  weights,
  r = 1,
  boot_type = c("bootstrap", "subbootstrap")
) {
  data[, paste0("bootrep_r_", r, "_", tag)] <- rep(0, nrow(data))
  data[which(data[, set] == tag), paste0("bootrep_r_", r, "_", tag)] <- 1

  nh0 <- table(data[, strata])
  if (boot_type == "bootstrap") {
    new_nh <- nh0
  } else {
    if (boot_type == "subbootstrap") {
      new_nh <- nh0 - 1
    }
  }

  nh0_tag <- table(data[which(data[,set] == tag), strata])

  for (hh in 1:length(unique(data[,strata]))) {
    if (nh0_tag[hh] < new_nh[hh]) {
      n_add  <- new_nh[hh] - nh0_tag[hh]
      id_opt <- which(data[, set] == tag & data[, strata] == hh)
      if (length(id_opt) > 1) {
        selected_id <- sample(id_opt, size = n_add, replace = TRUE)
      } else {
        selected_id <- rep(id_opt, n_add)
      }
      n_adds <- table(selected_id)
      data[as.numeric(names(table(selected_id))),
           paste0("bootrep_r_", r, "_", tag)] <- data[as.numeric(names(table(selected_id))), paste0("bootrep_r_", r, "_", tag)] + n_adds
    }
  }

  coef_h <- nh0/new_nh
  coef   <- coef_h[match(data[,strata], names(coef_h))]

  data[,paste0("rw_r_",r,"_", tag)] <- data[, weights] * data[, paste0("bootrep_r_", r, "_", tag)] * coef

  # Delete bootstrap repetition columns
  col_bootrep <- grep("bootrep_", colnames(data))
  data        <- data[,-col_bootrep]

  return(data)
}

# Extrapolation -----------------------------------------------------------

#' Extrapolation
#' @description This function allows calculating extrapolation.
#' @param data A data frame with information of cluster and strata indicators and sampling weights.
#' @param train_prob Proportion of PSUs set as training set.
#' @param strata Name of the column indicating strata
#' @param weights Name of the column indicating sampling weights
#' @param seed seed
#' @param R Number of times the sample is split
#' @return This function returns a new data frame with new columns, each of them indicating replicate weights for different subsets.
#' @export
split_strata <- function(
  data,
  train_prob,
  strata = NULL,
  weights,
  seed = 1,
  R = 1
) {

  if(is.null(strata)) { stop("Extrapolation method cannot be applied if strata are not defined") }

  set.seed(seed)
  seeds <- stats::runif(R) * 1000

  h <- unique(data[,strata])

  for (r in 1:R) {
    set.seed(seeds[r])

    number_h <- floor(length(h)*train_prob)
    train_h  <- sample(1:length(h), number_h)

    h_split        <- vector(length(h))
    names(h_split) <- h

    h_split[train_h]  <- "train"
    h_split[-train_h] <- "test"

    data[, paste0("set_", r)] <- as.vector(h_split[match(data[, strata], names(h_split))])

    data[, paste0("rw_r_", r, "_train")]         <- rep(0, nrow(data))
    id_train                                     <- which(data[, paste0("set_", r)] == "train")
    data[id_train, paste0("rw_r_", r, "_train")] <- data[id_train, weights]

    data[, paste0("rw_r_", r, "_test")]        <- rep(0, nrow(data))
    id_test                                    <- which(data[, paste0("set_", r)] == "test")
    data[id_test, paste0("rw_r_", r, "_test")] <- data[id_test, weights]

  }

  return(data)
}
