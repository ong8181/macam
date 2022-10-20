#' @title A wrapper function for MDR S-map
#' @description \code{s_map_mdr} performs MDR S-map introduced in Chang et al. (2021) Ecology Letters.
#' @param block Data.frame contains time series data.
#' @param effect_var Character. Column name that indicates effect variable.
#' @param E Numeric. Embedding dimensions that will be tested.
#' @param tp Numeric. Forecasting time ahead.
#' @param uic_method Character. Specify UIC method. `optimal` or `marginal`. See https://github.com/yutakaos/rUIC for detail.
#' @param only_uic Logical. If `TRUE`, only UIC results will be returned.
#' @param n_ssr Numeric. The total number of embeddings examined.
#' @param k Numeric. The number of embeddings used to calculate ensemble distance.
#' @param evaluated_by Character. Specify the criteria used to select the top k-th embeddings. Currently only `rho` is supported.
#' @param max_delay Numeric. The maximum number of delay to be used. This applies only to `effect_var`, and this may equal to `E-1`. Potential causal variables and their delays are determiend by UIC.
#' @param save_smap_coefficients Logical. If `TRUE`, S-map coefficients will be saved.
#' @param random_seed Numeric. Random seed.
#' @return A list that contains predictions, statistics, and S-map coefficients (if `save_smap_coefficients = TRUE`)
#' @noRd
#' @keywords internal
#' @details
#' \itemize{
#'  \item{Chang, C.-W., Miki, T., Ushio, M., Ke, P.-J., Lu, H.-P., Shiah, F.-K. & et al. (2021) Reconstructing large interaction networks from empirical time series data. Ecology Letters, 24, 2763– 2774. https://doi.org/10.1111/ele.13897}
#' }
s_map_mdr <- function(block, effect_var, E = 0:10, tp = -4:4,
                      uic_method = "optimal", only_uic = FALSE,
                      n_ssr = 10000, k = floor(sqrt(n_ssr)),
                      evaluate_by = "rho", max_delay = 2,
                      save_smap_coefficients = FALSE,
                      random_seed = 1234) {
  # Set random seed
  set.seed(random_seed)

  # Retrieve colnames
  x_names <- colnames(block)

  # ---------------------------------------------------- #
  # Check input arguments
  # ---------------------------------------------------- #
  if (!all(unique(x_names) == x_names)) stop("\"block\" should have unique column names.")
  # "block" should be data.frame
  if (!is.data.frame(block)) stop("\"block\" should be data.frame.")
  # "evaluated_by" should be "rho".
  if (evaluate_by != "rho") {
    stop("\"mae\" or \"rmse\" are not currently supported for \"evaluated_by\". Please use \"rho\" instead.")
  }
  if (!is.numeric(E) | !is.numeric(tp)) stop("\"E\" and \"tp\" should be numeric.")


  # ---------------------------------------------------- #
  # Identify causal relationship using rUIC
  # ---------------------------------------------------- #
  # Determining best embedding dimension (Univariate simplex)
  simp_x <- rUIC::simplex(block, lib_var = effect_var, E = E, tau = 1, tp = 1)
  Ex <- simp_x[which.min(simp_x$rmse),"E"]

  ## Prepare an object to store the output
  uic_res <- data.frame(matrix(NA, nrow = length(x_names[x_names != effect_var]), ncol = 9))
  colnames(uic_res) <- c("E","tau","tp","nn","n_lib","n_pred","rmse","te","pval")
  rownames(uic_res) <- x_names[x_names != effect_var]

  ## Perform UIC for all pairs
  for (y_i in x_names[x_names != effect_var]) {
    time_start <- proc.time()
    if (uic_method == "optimal") {
      # Testing the effect of "y_i" on "effect_var" using uic.optimal()
      uic_xy <- rUIC::uic.optimal(block, lib_var = effect_var, tar_var = y_i, E = E, tau = 1, tp = tp)
    } else if (uic_method == "marginal") {
      # Testing the effect of "y_i" on "effect_var" using uic.marginal()
      uic_xy <- rUIC::uic.marginal(block, lib_var = effect_var, tar_var = y_i, E = E, tau = 1, tp = tp)
    } else {
      stop("Specify valid uic_method: \"optimal\" or \"marginal\"")
    }
    uic_res[y_i, ] <- uic_xy[which.max(uic_xy$te),]
    ## Output message
    time_used <- (proc.time() - time_start)[3]
    message(sprintf("Effects from variable %s to %s tested: %.2f sec elapsed", y_i, effect_var, time_used))
  }
  # Delete the temporal object
  rm(uic_xy)

  # If a user wants to stop here, set only_uic = TRUE
  if (only_uic) return(uic_res)


  # ---------------------------------------------------- #
  # Construct multiview embedding
  # ---------------------------------------------------- #
  # Choose significant causal variables
  ## Exclude p > 0.05 & tp >= 1
  uic_sig <- uic_res[uic_res$pval <= 0.05 & uic_res$tp < 1, ]
  causal_var <- rownames(uic_sig)
  causal_var_tp <- paste0(causal_var, "_tp", uic_sig$tp)

  ## Make lagged block that includes all variables with a specific tp
  block_lagged <- rEDM::make_block(block[,effect_var], max_lag = max_delay + 1)[,-1]
  colnames(block_lagged) <- sprintf("%s_tp%d", effect_var, seq(0, - max_delay, by = -1))

  ## Generate lagged block of causal variables
  cause_lagged <- data.frame(matrix(NA, nrow = nrow(block_lagged), ncol = nrow(uic_sig)))
  for (i in 1:nrow(uic_sig)) {
    cause_i <- rownames(uic_sig)[i]
    tp_i <- uic_sig[i,"tp"]
    cause_lagged[,i] <- dplyr::lag(block[,cause_i], n = abs(tp_i))
    colnames(cause_lagged)[i] <- causal_var_tp[i]
  }
  block_multiview <- cbind(block_lagged, cause_lagged)
  valid_idx <- which(colnames(block_multiview) != paste0(effect_var, "_tp0"))

  ## Generate embedding list using sample() function
  potential_embeddings_list <- t(utils::combn(valid_idx, Ex-1, simplify = TRUE))

  if (n_ssr < nrow(potential_embeddings_list)) {
    embedding_idx <- sample(1:nrow(potential_embeddings_list), n_ssr)
    cause_var_embedding_list <- potential_embeddings_list[embedding_idx,]
  } else {
    cause_var_embedding_list <- potential_embeddings_list
  }


  # ---------------------------------------------------- #
  # Do multiview-embedded simplex projection using block_lnlp()
  # ---------------------------------------------------- #
  ## Prepare the result object of multiview analysis
  multiview_colnames <- c("embedding", "tp", "nn",
                          "num_pred", "rho", "mae",
                          "rmse", "perc", "p_val",
                          "const_pred_num_pred",
                          "const_pred_rho",
                          "const_pred_mae",
                          "const_pred_rmse",
                          "const_pred_perc",
                          "const_p_val")
  multiview_res <- data.frame(matrix(NA,
                                     ncol = length(multiview_colnames),
                                     nrow = nrow(cause_var_embedding_list)))
  colnames(multiview_res) <- multiview_colnames

  ## Main loop of randomized simplex projections
  for (i in 1:nrow(cause_var_embedding_list)) {
    ## Choose actual embedding ids (id = 1 is always a tarted column)
    multiview_idx <- c(1, cause_var_embedding_list[i,])
    ## Perform simplex projection
    multiview_res[i,] <- rEDM::block_lnlp(block_multiview[,multiview_idx],
                                          lib = c(1, nrow(block_multiview)),
                                          pred = c(1, nrow(block_multiview)),
                                          method = "simplex",
                                          silent = TRUE,
                                          tp = 1, target_column = 1)
    ## Replace embedding id
    multiview_res[i,]$embedding <- stringr::str_sub(stringr::str_c(multiview_idx, ",", collapse = " "), end = -2)
  }

  ## Sort the result based on "evaluated_by"
  multiview_res <- multiview_res[order(multiview_res[,evaluate_by], decreasing = T),]
  if (k < nrow(multiview_res)) {
    top_multiview_res <- multiview_res[1:k,]
  } else {
    top_multiview_res <- multiview_res
  }
  # Rename row names
  rownames(top_multiview_res) <- 1:nrow(top_multiview_res)


  # ---------------------------------------------------- #
  # Calculate ensemble weights for S-map
  # ---------------------------------------------------- #
  # Calculate a weight for each embedding
  w_multi <- top_multiview_res[,evaluate_by]/sum(top_multiview_res[,evaluate_by])

  # "multiview_dist" will be used to calculate S-map weights
  multiview_dist <- matrix(0, ncol = nrow(block), nrow = nrow(block))
  for (i in 1:nrow(top_multiview_res)) {
    top_multiview_ids <- stringr::str_split(top_multiview_res$embedding[i], pattern = ",")[[1]] %>%
      as.numeric
    multiview_dist <- multiview_dist + w_multi[i] * as.matrix(stats::dist(block_multiview[,top_multiview_ids]))
  }


  # ---------------------------------------------------- #
  # Perform MDR S-map using macam::extended_lnlp
  # ---------------------------------------------------- #
  mdr_res <- macam::extended_lnlp(block_multiview,
                                  lib = 1:nrow(block_multiview),
                                  pred = 1:nrow(block_multiview),
                                  target_column = 1,
                                  dist_w = multiview_dist,
                                  theta = 1,
                                  regularized = TRUE,
                                  lambda = 0.1,
                                  alpha = 0,
                                  save_smap_coefficients = save_smap_coefficients)

  # Return results
  return(mdr_res)
}
