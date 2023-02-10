#' @title Perform UIC across data.frame
#' @description \code{uic_across} Perform UIC for a target variable and multiple causal variables
#' @param block Data.frame contains time series data.
#' @param effect_var Character or Numeric. Column name or index of the effect variable.
#' @param uic_method Character. `optimal` or `marginal`. For detail, see https://github.com/yutakaos/rUIC.
#' @param E_range Numeric. Embedding dimensions that will be tested.
#' @param tp_range Numeric. `tp` tested for UIC.
#' @param random_seed Numeric. Random seed.
#' @param silent Logical. if `TRUE`, progress message will not be shown.
#' @return A data.frame that contains UIC results
#' @details
#' \itemize{
#'  \item{Osada & Ushio (2020) rUIC:Unified Information-theoretic Causality for R. https://doi.org/10.5281/zenodo.5163234}
#' }
#' @export
uic_across <- function(block,
                       effect_var,
                       uic_method = "optimal",
                       E_range = 0:10,
                       tp_range = -4:0,
                       random_seed = 1234,
                       silent = FALSE) {
  # Set random seed
  #set.seed(random_seed)

  # Retrieve colnames
  x_names <- colnames(block)

  # Check input arguments
  if (!all(unique(x_names) == x_names)) stop("\"block\" should have unique column names.")
  if (!is.data.frame(block)) stop("\"block\" should be data.frame.")
  if (!is.numeric(E_range) | !is.numeric(tp_range)) stop("\"E_range\" and \"tp_range\" should be numeric.")
  if (uic_method != "optimal" & uic_method != "marginal") stop("\"uic_method\" should be \"optimal\" or \"marginal\".")
  if (is.numeric(effect_var)) effect_var <- x_names[effect_var]

  # ---------------------------------------------------- #
  # Preparation for UIC
  # ---------------------------------------------------- #
  # Rearrange the column order
  block <- dplyr::select(block, tidyselect::all_of(effect_var), dplyr::everything())

  # ---------------------------------------------------- #
  # Identify causal relationship using rUIC
  # ---------------------------------------------------- #
  ## Perform UIC for all pairs
  if (uic_method == "optimal") {
    # ---------------------------------------------------- #
    # uic.optimal()
    # ---------------------------------------------------- #
    for (y_i in x_names[x_names != effect_var]) {
      time_start <- proc.time()
      # Testing the effect of "y_i" on "effect_var" using uic.optimal()
      uic_xy <- rUIC::uic.optimal(block, lib_var = effect_var, tar_var = y_i, E = E_range, tau = 1, tp = tp_range, seed = random_seed) %>%
        dplyr::mutate(effect_var = effect_var, cause_var = y_i)
      # Combine results
      if (y_i != x_names[x_names != effect_var][1]) { uic_res <- rbind(uic_res, uic_xy) } else { uic_res <- uic_xy }
      # Output message
      time_used <- (proc.time() - time_start)[3]
      if (!silent) { message(sprintf("Effects from %s to %s tested by UIC: %.2f sec elapsed", y_i, effect_var, time_used)) }
    }
  } else if (uic_method == "marginal") {
    # ---------------------------------------------------- #
    # uic.marginal()
    # ---------------------------------------------------- #
    for (y_i in x_names[x_names != effect_var]) {
      time_start <- proc.time()
      # Testing the effect of "y_i" on "effect_var" using uic.marginal()
      uic_xy <- rUIC::uic.marginal(block, lib_var = effect_var, tar_var = y_i, E = E_range, tau = 1, tp = tp_range, seed = random_seed) %>%
        dplyr::mutate(effect_var = effect_var, cause_var = y_i)
      # Combine results
      if (y_i != x_names[x_names != effect_var][1]) { uic_res <- rbind(uic_res, uic_xy) } else { uic_res <- uic_xy }
      # Output message
      time_used <- (proc.time() - time_start)[3]
      if (!silent) { message(sprintf("Effects from %s to %s tested by UIC: %.2f sec elapsed", y_i, effect_var, time_used)) }
    }
  }

  # Message
  if (!silent) {
    message("By using `uic_across()` the effect variable is recorded in the column named `effect_var`, while potential causal variables are recorded in the column named `cause_var`.")
    message("Please note that these column names are used in the subsequent analysis if you use `make_block_mvd()`.")
  }

  # Return results
  return(as.data.frame(uic_res))
}



#' @title Generate block for MDR S-map
#' @description \code{make_block_mvd} Generate data.frame for MDR S-map
#' @param block Data.frame contains time series data.
#' @param uic_res Data.frame contains UIC results. Usually this is the output of `uic_across()`, but custom data.frame may be usable.
#' @param effect_var Character or Numeric. Column name or index of the effect variable.
#' @param E_effect_var Numeric. Optimal embedding dimension of the effect variable.
#' @param cause_var_colname Character. A column name for causal variables in `uic_res`.
#' @param include_var Character. `all_significant`, `strongest_only`, or `tp0_only`. If `all_significant`, all significantly influencing variables are used for the embedding. If `strongest_only`, tp with the strongest influence for each variable is used. If `tp0_only`, only variables with no time-delay are used.
#' @param p_threshold Numeric. Rows with pval larger than `p_threshold` is removed.
#' @param sort_tp Logical. If TRUE, `block_mvd` is sorted according to `tp` for each causal variable.
#' @param silent Logical. if `TRUE`, progress message will not be shown.
#' @return Embedded time series will be returned.
#' @export
make_block_mvd <- function (block,
                            uic_res,
                            effect_var,
                            E_effect_var,
                            cause_var_colname = "cause_var",
                            include_var = "strongest_only",
                            p_threshold = 0.050,
                            sort_tp = TRUE,
                            silent = FALSE) {
  # Retrieve colnames
  x_names <- colnames(block)

  # Check colnames to be used (and critical) in the analysis
  if (is.numeric(effect_var)) effect_var <- x_names[effect_var]
  if (!(effect_var %in% x_names)) stop("No 'effect_var' in the block!")
  if (!(cause_var_colname %in% colnames(uic_res))) stop("No 'cause_var_colname' in `uic_res`! Please specify the correct colname for causal variables")
  if (!("tp" %in% colnames(uic_res))) stop("'tp' column is required for `uic_res`.")
  if (!("te" %in% colnames(uic_res))) stop("'te' column is required for `uic_res`.")
  if (!("pval" %in% colnames(uic_res))) stop("'pval' column is required for `uic_res`.")

  # Check arguments
  if (!all(unique(x_names) == x_names)) stop("\"block\" should have unique column names.")
  if (!is.data.frame(block)) stop("\"block\" should be data.frame.")
  if (!is.numeric(p_threshold)) stop("\"p_threshold\" should be numeric.")
  if (!include_var %in% c("all_significant", "strongest_only", "tp0_only")) stop("\"include_var\" should be \"all_significant\", \"strongest_only\", or \"tp0_only\".")
  if (E_effect_var < 1) stop("\"max_delay_self\" should be >= 1.")

  # Preparation
  if (include_var == "tp0_only") {
    block_mvd <- data.frame(block[,effect_var])
    colnames(block_mvd) <- sprintf("%s_tp0", effect_var)
  } else {
    block_mvd <- data.frame(rEDM::make_block(block[,effect_var], max_lag = E_effect_var)[,-1])
    colnames(block_mvd) <- sprintf("%s_tp%s", effect_var, 0:(-(E_effect_var-1)))
  }
  # Pre-screening (p & tp)
  if (!silent) message("UIC results with `tp` <= 0 and `pval` <= 0.05 are kept for further analyses.")
  uic_res <- uic_res[uic_res$pval <= p_threshold & uic_res$tp <= 0,]
  if (nrow(uic_res) < 1) stop("No significant causal variables were detected. Please use the univariate S-map.")

  # Sort block columns according to tp for each causal variable
  if (sort_tp) {
    sort_id <- 0
    for (col_i in unique(uic_res[,cause_var_colname])) {
      sort_id_i <- order(-uic_res[uic_res[,cause_var_colname] == col_i,"tp"])
      sort_id <- c(sort_id, max(sort_id) + sort_id_i)
    }
    sort_id <- sort_id[-1]; uic_res <- uic_res[sort_id,]
  }

  # ---------------------------------------------------- #
  # Select variables to be included in block
  # ---------------------------------------------------- #
  if (include_var == "all_significant") {
    # ---------------------------------------------------- #
    # Select all significant variables
    # ---------------------------------------------------- #
    for (i in 1:nrow(uic_res)) {
      block_new <- dplyr::lag(block[,uic_res[i, cause_var_colname]], n = abs(uic_res[i,"tp"]))
      block_new <- data.frame(block_new)
      colnames(block_new) <- sprintf("%s_tp%s", uic_res[i, cause_var_colname], uic_res[i,"tp"])
      block_mvd <- cbind(block_mvd, block_new)
    }
  } else if (include_var == "strongest_only") {
    # ---------------------------------------------------- #
    # Select tp with the strongest influence from each causal variable
    # (te is used as a criterion)
    # ---------------------------------------------------- #
    for (cause_i in unique(uic_res[,cause_var_colname])) {
      uic_res_tmp <- uic_res[uic_res[,cause_var_colname] == cause_i,]
      if (!exists("uic_res_new")) {
        uic_res_new <- uic_res_tmp[which.max(uic_res_tmp[,"te"]),]
      } else {
        uic_res_new <- rbind(uic_res_new, uic_res_tmp[which.max(uic_res_tmp[,"te"]),])
      }
    }
    # Replace "uic_res"
    uic_res <- uic_res_new
    #uic_res <- uic_res %>% dplyr::group_by(.data$cause_var) %>% dplyr::filter(.data$te == max(.data$te))
    for (i in 1:nrow(uic_res)) {
      block_new <- dplyr::lag(block[,uic_res[i,cause_var_colname]], n = abs(uic_res[i,"tp"]))
      block_new <- data.frame(block_new)
      colnames(block_new) <- sprintf("%s_tp%s", uic_res[i,cause_var_colname], uic_res[i,"tp"])
      block_mvd <- cbind(block_mvd, block_new)
    }
  } else if (include_var == "tp0_only") {
    # ---------------------------------------------------- #
    # Select causal variables with tp = 0
    # ---------------------------------------------------- #
    block_new <- block[,unique(uic_res[,cause_var_colname])]
    block_new <- data.frame(block_new)
    colnames(block_new) <- sprintf("%s_tp0", unique(uic_res[,cause_var_colname]))
    block_mvd <- cbind(block_mvd, block_new)
  }
  # Return results
  return(block_mvd)
}



#' @title Computing multiview distance
#' @description \code{compute_mvd} Compute multiview distance
#' @param block_mvd Data.frame contains time series data. The first column should be the target column.
#' @param effect_var Character or Numeric. Column name or index of the effect variable.
#' @param E Numeric. Optimal embedding dimension of `effect_var`
#' @param tp Numeric. Forecasting time ahead.
#' @param lib Numeric vector. Library indices.
#' @param pred Numeric vector. Prediction indices.
#' @param make_block_method Character. If `naive`, the multivariate data.frame (`block_mvd`) is directly used to calculate the multiview distance. If `rEDM`, `rEDM::make_block()` is used to add the time-delayed ordinate for each variable, which make the method equivalent to Chang et al. (2021).
#' @param make_block_max_lag Numeric. `max_lag` in `rEDM::make_block()`. This argument will be used only if `make_block_method = "rEDM"`.
#' @param n_ssr Numeric. The total number of embeddings examined.
#' @param k Numeric. The number of embeddings used to calculate ensemble distance.
#' @param random_seed Numeric. Random seed.
#' @param distance_only Logical. if `TRUE`, only distance matrix is returned.
#' @param silent Logical. if `TRUE`, progress message will not be shown.
#' @return A distance matrix and other information (if `distance_only = FALSE`).
#' @details
#' \itemize{
#'  \item{Chang et al. (2021) Ecology Letters. https://doi.org/10.1111/ele.13897}
#' }
#' @export
compute_mvd <- function (block_mvd, effect_var, E,
                         tp = 1,
                         lib = c(1, nrow(block_mvd)),
                         pred = lib,
                         make_block_method = "naive", # "rEDM"
                         make_block_max_lag = E,
                         n_ssr = 10000, k = floor(sqrt(n_ssr)),
                         random_seed = 1234,
                         distance_only = TRUE,
                         silent = FALSE) {
  # Set random seed
  set.seed(random_seed)
  # Retrieve colnames
  x_names <- colnames(block_mvd)

  # Check input
  if (!(E %% 1 == 0) | E < 1) stop("\"E\" should be a natural number.")
  if (is.numeric(effect_var)) effect_var <- x_names[effect_var]
  if (colnames(block_mvd)[1] != sprintf("%s_tp0", effect_var)) stop("The first column should be the target variable with tp = 0!")
  if (!(make_block_method %in% c("naive", "rEDM"))) stop("\"make_block_method\" should be \"naive\" or \"rEDM\".")


  # ---------------------------------------------------- #
  # Generate an embedding list
  # ---------------------------------------------------- #
  if (E == 1) {
    embedding_list <- matrix(1, ncol = 1, nrow = 1)
  } else {
    if (make_block_method == "naive") {
      # ---------------------------------------------------- #
      # Generate embedding list using sample() function
      # ---------------------------------------------------- #
      valid_idx <- 2:length(x_names)
    } else if (make_block_method == "rEDM") {
      # ---------------------------------------------------- #
      # ! Recommended only if "block_mvd" is generated by specifying "tp0_only" in make_block_mvd()
      # Make a new block using rEDM::make_block()
      # ---------------------------------------------------- #
      block_redm <- rEDM::make_block(block_mvd, max_lag = make_block_max_lag)[,-1]
      block_mvd <- block_redm
      valid_idx <- 2:ncol(block_mvd)
    }
    # Check the number of possible embeddings
    if (choose(length(valid_idx), E-1) >= n_ssr) {
      # Output message
      if (!silent) message(sprintf("The performance of %s embeddings (n_ssr) will be evaluated from %s potential embeddings", n_ssr, choose(length(valid_idx), E-1)))
      # Randomly select embeddings
      emb_list0 <- matrix(0, nrow = 1, ncol = E-1); i <- 1
      while (nrow(emb_list0) < n_ssr + 1) {
        i <- i + 1; set.seed(random_seed + i)
        emb_list0_tmp <- matrix(sort(sample(valid_idx, E-1, replace = FALSE)), nrow = 1)
        emb_dif <- t(t(emb_list0) - c(emb_list0_tmp)) # Subtract a vector from matrix
        if(all(!apply(emb_dif == 0, 1, all))) emb_list0 <- rbind(emb_list0, emb_list0_tmp)
      }
      # Add target variable
      emb_list1 <- cbind(matrix(1, ncol = 1, nrow = nrow(emb_list0)-1), emb_list0[-1,])
      # Sort emb_list1
      emb_list1 <- emb_list1[do.call(order, as.list(as.data.frame(emb_list1))),]
      embedding_list <- emb_list1
    } else {
      if (!silent) message("The number of potential embeddings is smaller than \"n_srr\". Thus, the performance of all potential embeddings will be evalulated.")
      emb_list0 <- t(utils::combn(valid_idx, E-1, simplify = TRUE))
      emb_list1 <- cbind(matrix(1, ncol = 1, nrow = nrow(emb_list0)), emb_list0)
      embedding_list <- emb_list1
    }
  }
  nrow_res <- nrow(embedding_list)

  # ---------------------------------------------------- #
  # Main loop of randomized simplex projections
  # ---------------------------------------------------- #
  for (i in 1:nrow_res) {
    ## Choose embedding ids (id = 1 is always a target column)
    embedding_idx_i <- embedding_list[i,]
    ## Perform simplex projection
    rand_embed_res_i <- rEDM::block_lnlp(block_mvd[,embedding_idx_i],
                                         lib = lib,
                                         pred = pred,
                                         method = "simplex",
                                         silent = TRUE,
                                         tp = tp,
                                         target_column = 1)
    # Add information
    rand_embed_res_i$embedding <- stringr::str_sub(stringr::str_c(embedding_idx_i, ",", collapse = " "), end = -2)
    # Merge results
    if (i == 1) { rand_embed_res <- rand_embed_res_i } else { rand_embed_res <- rbind(rand_embed_res, rand_embed_res_i)}
  }

  ## Remove rho <= 0 because it is meaningless prediction.
  rand_embed_res <- rand_embed_res[rand_embed_res[,"rho"] > 0,]
  if (!silent) message("Embeddings with rho < 0 were removed.")
  if (nrow(rand_embed_res) == 0) stop("Cannot find embeddings with positive prediction skill. Reconsider parameters.")

  ## Sort the result based on "evaluated_by"
  rand_embed_res <- rand_embed_res[order(rand_embed_res[,"rho"], decreasing = T),]
  if (k <= nrow(rand_embed_res)) {
    # Output message
    if (!silent) message(sprintf("Top %s embeddings (k) will be used from %s embeddings (rho > 0 in n_ssr embeddings)", k, nrow(rand_embed_res)))
    top_embed_res <- rand_embed_res[1:k,]
  } else {
    if (!silent) message("The total number of valid embeddings is smaller than \"k\". Thus, all valid embeddings will be used to calculate the multiview distance.")
    top_embed_res <- rand_embed_res
  }
  # Rename row names
  rownames(top_embed_res) <- 1:nrow(top_embed_res)

  # ---------------------------------------------------- #
  # Calculate ensemble weights for S-map
  # ---------------------------------------------------- #
  # Calculate a weight for each embedding
  w_multi <- top_embed_res[,"rho"]/sum(top_embed_res[,"rho"])

  # "multiview_dist" will be used to calculate S-map weights
  multiview_dist <- matrix(0, ncol = nrow(block_mvd), nrow = nrow(block_mvd))
  for (i in 1:nrow(top_embed_res)) {
    top_embed_ids <- stringr::str_split(top_embed_res$embedding[i], pattern = ",")[[1]] %>%
      as.numeric
    multiview_dist <- multiview_dist + w_multi[i] * as.matrix(stats::dist(block_mvd[,top_embed_ids]))
  }

  # ---------------------------------------------------- #
  # Return results
  # ---------------------------------------------------- #
  if (distance_only) {
    return(multiview_dist)
  } else {
    all_res <- list(multiview_dist); names(all_res) <- "multiview_dist"
    all_res$embeddings <- rand_embed_res
    all_res$top_embeddings <- top_embed_res
    if (is.null(random_seed)) random_seed <- "NULL"
    all_res$parms <- data.frame(E = E, tp = tp, n_ssr = nrow(rand_embed_res),
                                k = nrow(top_embed_res), random_seed = random_seed)
    return(all_res)
  }
}


#' @title A function for MDR S-map
#' @description \code{s_map_mdr} performs MDR S-map introduced in Chang et al. (2021) Ecology Letters.
#' @param block_mvd Data.frame contains time series data. The first column should be the target column.
#' @param dist_w Matrix. Multiview disntance matrix. Usually this is the output of `compute_mvd()`.
#' @param lib Numeric vector. Library indices.
#' @param pred Numeric vector. Prediction indices.
#' @param tp Numeric. Forecasting time ahead.
#' @param theta Numeric. Weighing function for S-map.
#' @param weight_method Specify weighing method for `dist_w`. Default is `sqrt`, which is the same for Chang et al. (2021).
#' @param regularized Logical If `TRUE`, regularized S-map will be performed. If `FALSE`, the normal S-map will be performed. Please use `rEDM::s_map` function.
#' @param lambda Numeric. Specify the strength of penalty in the regularization.
#' @param alpha Numeric. `alpha = 0` is the ridge regression, `alpha = 1` is the lasso regression, and `0 < alpha < 1` is an elastic net.
#' @param glmnet_parallel Logical. If TRUE, the computation will be parallel (currently, experimental).
#' @param save_smap_coefficients Logical. If `TRUE`, S-map coefficients will be saved.
#' @param random_seed Numeric. Random seed.
#' @return A list that contains predictions, statistics, and S-map coefficients (if `save_smap_coefficients = TRUE`)
#' @details
#' \itemize{
#'  \item{Chang, C.-W., Miki, T., Ushio, M., Ke, P.-J., Lu, H.-P., Shiah, F.-K. & et al. (2021) Reconstructing large interaction networks from empirical time series data. Ecology Letters, 24, 2763– 2774. https://doi.org/10.1111/ele.13897}
#' }
#' @export
s_map_mdr <- function(block_mvd,
                      dist_w,
                      lib = c(1, nrow(block_mvd)), pred = lib,
                      tp = 1, theta = 8,
                      weight_method = "sqrt",
                      regularized = FALSE,
                      lambda = 0, alpha = 0,
                      glmnet_parallel = FALSE,
                      save_smap_coefficients = FALSE,
                      random_seed = 1234) {
  # Check input
  if (!(nrow(block_mvd) == nrow(dist_w))) stop(" \"dist_w\" should have the same size with \"block_mvd\"")
  if (!is.null(weight_method) & weight_method != "sqrt") stop("Specify the valid option for \"weight_method\": NULL or \"sqrt\"")
  if (weight_method == "sqrt") dist_w = sqrt(dist_w)

  # ---------------------------------------------------- #
  # Perform MDR S-map using macam::extended_lnlp
  # ---------------------------------------------------- #
  mdr_res <- macam::extended_lnlp(block_mvd,
                                  lib = lib,
                                  pred = pred,
                                  target_column = 1,
                                  tp = tp,
                                  dist_w = dist_w,
                                  theta = theta,
                                  regularized = regularized,
                                  lambda = lambda,
                                  alpha = alpha,
                                  glmnet_parallel = glmnet_parallel,
                                  save_smap_coefficients = save_smap_coefficients,
                                  random_seed = random_seed)
  # Return results
  return(mdr_res)
}



#' @title A all-in-one wrapper function for MDR S-map
#' @description \code{s_map_mdr_all} A all-in-one wrapper function for MDR S-map. Detect causality between variables, construct data.frame for multiview embedding, compute multiview distance, and perform MDR S-map. For fine-tuning, use the step-by-step functions such as `compute_mvd()` etc..
#' @param block Data.frame contains time series data. The first column should be the target column.
#' @param effect_var Character or Numeric. Column name or index of the effect variable.
#' @param lib Numeric vector. Library indices.
#' @param pred Numeric vector. Prediction indices.
#' @param tp Numeric. Forecasting time ahead.
#' @param E_range Numeric. Embedding dimensions that will be tested.
#' @param tp_range Numeric. `tp` tested for UIC.
#' @param n_ssr Numeric. The total number of embeddings examined.
#' @param k Numeric. The number of embeddings used to calculate ensemble distance.
#' @param theta Numeric. Weighing function for S-map.
#' @param regularized Logical If `TRUE`, regularized S-map will be performed. If `FALSE`, the normal S-map will be performed. Please use `rEDM::s_map` function.
#' @param lambda Numeric. Specify the strength of penalty in the regularization.
#' @param alpha Numeric. `alpha = 0` is the ridge regression, `alpha = 1` is the lasso regression, and `0 < alpha < 1` is an elastic net.
#' @param glmnet_parallel Logical. If TRUE, the computation will be parallel (currently, experimental).
#' @param save_smap_coefficients Logical. If `TRUE`, S-map coefficients will be saved.
#' @param silent Logical. if `TRUE`, progress message will not be shown.
#' @param random_seed Numeric. Random seed.
#' @return A list that contains predictions, statistics, and S-map coefficients (if `save_smap_coefficients = TRUE`)
#' @export
s_map_mdr_all <- function (block,
                           effect_var,
                           lib = c(1, nrow(block)),
                           pred = lib,
                           tp = 1,
                           E_range = 0:10,
                           tp_range = -4:0,
                           n_ssr = 100,
                           k = 10,
                           theta = 0,
                           regularized = FALSE,
                           lambda = 0,
                           alpha = 0,
                           glmnet_parallel = FALSE,
                           save_smap_coefficients = TRUE,
                           silent = FALSE,
                           random_seed = 1234) {
  # Step 1: Determine optimal embedding dimension
  simp_res <- rUIC::simplex(block, lib_var = effect_var, E = E_range, tau = 1, tp = 1)
  Ex <- simp_res[which.min(simp_res$rmse),"E"]

  # Step 2: Perform UIC to detect causality
  uic_res <- uic_across(block, effect_var, E_range = E_range, tp_range = tp_range, silent = silent, random_seed = random_seed)

  # Step 3: Make block to calculate multiview distance
  block_mvd <- make_block_mvd(block, uic_res, effect_var, E_effect_var = Ex, include_var = "tp0_only")

  # Step. 4: Compute multiview distance
  multiview_dist <- compute_mvd(block_mvd, effect_var, E = Ex,
                                make_block_method = "rEDM",
                                make_block_max_lag = Ex,
                                n_ssr = n_ssr, k = k,
                                tp = tp, random_seed = random_seed)

  # Step. 5: Do MDR S-map
  mdr_res <- s_map_mdr(block_mvd,
                       dist_w = multiview_dist,
                       lib = lib,
                       pred = pred,
                       tp = tp,
                       theta = theta,
                       regularized = regularized,
                       lambda = lambda,
                       alpha = alpha,
                       glmnet_parallel = glmnet_parallel,
                       save_smap_coefficients = save_smap_coefficients,
                       random_seed = random_seed)

  return(mdr_res)
}


