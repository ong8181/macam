#' @title A wrapper function for MDR S-map
#' @description \code{s_map_mdr} performs MDR S-map introduced in Chang et al. (2021) Ecology Letters.
#' @param block Data.frame contains time series data.
#' @param effect_var Character or Numeric. Column name or index of the effect variable.
#' @param lib Numeric vector. Library indices.
#' @param pred Numeric vector. Prediction indices.
#' @param E_var Numeric. If specified, `E_var` is used as an optimal embedding dimension of the effect variable. If `NULL`, `rUIC::simplex()` and `E_range` is used to determine the optimal embedding dimension. Default is `NULL`.
#' @param E_range Numeric. Embedding dimensions that will be tested.
#' @param tp Numeric. Forecasting time ahead.
#' @param tp_range Numeric. `tp` tested for UIC.
#' @param uic_method Character. Specify UIC method. `optimal` or `marginal`. See https://github.com/yutakaos/rUIC for detail.
#' @param uic_res Data.frame. Should contain UIC results. Required if `anaysis_flow = "mdr"`. Colnames of `uic_res` should be `c("E","tau","tp","nn","n_lib","n_pred","rmse","te","pval")`.
#' @param analysis_flow Character or numeric. `uic` or `1` performs UIC analysis only, and `mdr` or `2` performs MDR S-map only.  `uic+mdr` or `3` performs UIC + MDR S-map. Default is `uic+mdr`.
#' @param n_ssr Numeric. The total number of embeddings examined.
#' @param k Numeric. The number of embeddings used to calculate ensemble distance.
#' @param evaluate_by Character. Specify the criteria used to select the top k-th embeddings. Currently only `rho` is supported.
#' @param use_all_lag Logical. If `TRUE`, use all lagged variables (< max_lag or E_var) of the causal variable. For example, if a variable `A` at time lag `-2` is causal and if `max_lag = 3`, then `A_tp-2`, `A_tp-3`, `A_tp-4`, and `A_tp-5` will be used for the multiview embedding.
#' @param max_delay Numeric. The maximum number of delay to be used. Default is the optimal embedding dimension (E) of `effect_var` minus 1. This applies only to `effect_var`. E of potential causal variables and their delays are determiend by UIC.
#' @param theta Numeric. Weighing function for S-map.
#' @param regularized Logical If `TRUE`, regularized S-map will be performed. If `FALSE`, the normal S-map will be performed. Please use `rEDM::s_map` function.
#' @param lambda Numeric. Specify the strength of penalty in the regularization.
#' @param alpha Numeric. `alpha = 0` is the ridge regression, `alpha = 1` is the lasso regression, and `0 < alpha < 1` is an elastic net.
#' @param glmnet_parallel Logical. If TRUE, the computation will be parallel (currently, experimental).
#' @param save_smap_coefficients Logical. If `TRUE`, S-map coefficients will be saved.
#' @param random_seed Numeric. Random seed.
#' @param silent Logical. if `TRUE`, progress message will not be shown.
#' @return A list that contains predictions, statistics, and S-map coefficients (if `save_smap_coefficients = TRUE`)
#' @details
#' \itemize{
#'  \item{Chang, C.-W., Miki, T., Ushio, M., Ke, P.-J., Lu, H.-P., Shiah, F.-K. & et al. (2021) Reconstructing large interaction networks from empirical time series data. Ecology Letters, 24, 2763– 2774. https://doi.org/10.1111/ele.13897}
#' }
#' @export
s_map_mdr <- function(block,
                      effect_var,
                      lib = c(1, nrow(block)),
                      pred = lib,
                      E_var = NULL,
                      E_range = 0:10,
                      tp = 1,
                      tp_range = -4:0,
                      uic_method = "optimal",
                      uic_res = NULL,
                      analysis_flow = "uic+mdr",
                      n_ssr = 10000,
                      k = floor(sqrt(n_ssr)),
                      evaluate_by = "rho",
                      use_all_lag = FALSE,
                      max_delay = NULL,
                      theta = 0,
                      regularized = FALSE,
                      lambda = NULL,
                      alpha = 0, # default is the ridge regression.
                      # If alpha = 1, it is the lasso regression
                      glmnet_parallel = FALSE,
                      save_smap_coefficients = FALSE,
                      random_seed = 1234,
                      silent = FALSE) {
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
  if (!is.numeric(E_range) | !is.numeric(tp_range)) stop("\"E_range\" and \"tp_range\" should be numeric.")
  if (is.numeric(effect_var)) effect_var <- x_names[effect_var]

  # ---------------------------------------------------- #
  # Choose analysis flow
  # ---------------------------------------------------- #
  if (analysis_flow == "uic+mdr" | analysis_flow == "uic +mdr" | analysis_flow == "uic+ mdr" |
      analysis_flow == "uic + mdr" | analysis_flow == "all" | analysis_flow == "both" | analysis_flow == 3) {
    analysis_flow <- "uic+mdr"
  } else if (analysis_flow == "uic" | analysis_flow == "uic only" | analysis_flow == "uic_only" | analysis_flow == 1) {
    analysis_flow <- "uic"
  } else if (analysis_flow == "mdr" | analysis_flow == "mdr only" | analysis_flow == "mdr_only" | analysis_flow == 2) {
    analysis_flow <- "mdr"
  } else {
    stop("Specify valid \"analysis_flow\": uic+mdr, uic, or mdr.")
  }


  # ---------------------------------------------------- #
  # Identify causal relationship using rUIC
  # ---------------------------------------------------- #
  # Rearrange the column order
  block <- dplyr::select(block, tidyselect::all_of(effect_var), dplyr::everything())

  # Determining best embedding dimension (Univariate simplex)
  if (is.null(E_var)) {
    simp_x <- rUIC::simplex(block, lib_var = effect_var, E = E_range, tau = 1, tp = 1)
    Ex <- simp_x[which.min(simp_x$rmse),"E"]
  } else {
    Ex <- E_var
  }
  if (!silent) message(sprintf("Note: E of the effect variable = %s", Ex))

  if (analysis_flow == "uic+mdr" | analysis_flow == "uic") {
    ## Prepare an object to store the output
    uic_res <- data.frame(matrix(NA, nrow = length(x_names[x_names != effect_var]), ncol = 9))
    colnames(uic_res) <- c("E","tau","tp","nn","n_lib","n_pred","rmse","te","pval")
    rownames(uic_res) <- x_names[x_names != effect_var]

    ## Perform UIC for all pairs
    for (y_i in x_names[x_names != effect_var]) {
      time_start <- proc.time()
      if (uic_method == "optimal") {
        # Testing the effect of "y_i" on "effect_var" using uic.optimal()
        uic_xy <- rUIC::uic.optimal(block, lib_var = effect_var, tar_var = y_i, E = E_range, tau = 1, tp = tp_range)
      } else if (uic_method == "marginal") {
        # Testing the effect of "y_i" on "effect_var" using uic.marginal()
        uic_xy <- rUIC::uic.marginal(block, lib_var = effect_var, tar_var = y_i, E = E_range, tau = 1, tp = tp_range)
      } else {
        stop("Specify valid uic_method: \"optimal\" or \"marginal\"")
      }
      uic_res[y_i, ] <- uic_xy[which.max(uic_xy$te),]
      ## Output message
      time_used <- (proc.time() - time_start)[3]
      if (!silent) {
        message(sprintf("Effects from %s to %s tested by UIC: %.2f sec elapsed", y_i, effect_var, time_used))
      }
    }
    # Delete the temporal object
    rm(uic_xy); rm(time_start); rm(time_used)

    # If a user wants to stop here, set analysis_flow == "uic" or "uic_only"
    if (analysis_flow == "uic") return(uic_res)
    }


  # ---------------------------------------------------- #
  # Construct multiview embedding
  # ---------------------------------------------------- #
  if (analysis_flow == "uic+mdr" | analysis_flow == "mdr") {
    # Choose significant causal variables
    ## Exclude p > 0.05 & tp >= 1
    uic_sig <- uic_res[uic_res$pval <= 0.05 & uic_res$tp < 1, ]
    if (nrow(uic_sig) < 1) {
      stop("No causal variables are detected. Perform univariate analysis instead.")
    }
    causal_var <- rownames(uic_sig)
    causal_var_tp <- paste0(causal_var, "_tp", uic_sig$tp)

    ## Make lagged block that includes all variables with a specific tp
    if (!is.null(max_delay)) {
      block_lagged <- data.frame(rEDM::make_block(block[,effect_var], max_lag = max_delay + 1)[,-1])
      colnames(block_lagged) <- sprintf("%s_tp%d", effect_var, seq(0, - max_delay, by = -1))
    } else {
      block_lagged <- data.frame(rEDM::make_block(block[,effect_var], max_lag = Ex)[,-1])
      colnames(block_lagged) <- sprintf("%s_tp%d", effect_var, seq(0, - (Ex-1), by = -1))
    }

    ## Generate lagged block of causal variables
    if (use_all_lag) {
      if (!is.null(max_delay)) {max_lag2 <- max_delay + 1} else {max_lag2 <- Ex}
      cause_lagged <- data.frame(matrix(NA, nrow = nrow(block_lagged),
                                        ncol = nrow(uic_sig)*max_lag2))
      for (i in 1:nrow(uic_sig)) {
        col_range <- ((i-1)*max_lag2+1):(i*max_lag2)
        cause_i <- rownames(uic_sig)[i]
        tp_i <- uic_sig[i,"tp"]
        # Make lagged block (Note: lagged from the causal tp. Might be modified in future)
        cause_lagged[,col_range] <- rEDM::make_block(dplyr::lag(block[,cause_i], n = abs(tp_i)) %>% data.frame,
                                                     max_lag = max_lag2)[,2:(max_lag2+1)]
        new_colnames_i <- sprintf("%s_tp%d", cause_i, tp_i - 0:(max_lag2-1))
        colnames(cause_lagged)[col_range] <- new_colnames_i
      }
      block_multiview <- cbind(block_lagged, cause_lagged)
      valid_idx <- which(colnames(block_multiview) != paste0(effect_var, "_tp0"))
    } else {
      cause_lagged <- data.frame(matrix(NA, nrow = nrow(block_lagged), ncol = nrow(uic_sig)))
      for (i in 1:nrow(uic_sig)) {
        cause_i <- rownames(uic_sig)[i]
        tp_i <- uic_sig[i,"tp"]
        cause_lagged[,i] <- dplyr::lag(block[,cause_i], n = abs(tp_i))
        colnames(cause_lagged)[i] <- causal_var_tp[i]
      }
      block_multiview <- cbind(block_lagged, cause_lagged)
      valid_idx <- which(colnames(block_multiview) != paste0(effect_var, "_tp0"))
    }

    ## Generate embedding list using sample() function
    if (Ex > 1) {
      potential_embeddings_list <- t(utils::combn(valid_idx, Ex-1, simplify = TRUE))

      if (n_ssr < nrow(potential_embeddings_list)) {
        embedding_idx <- sample(1:nrow(potential_embeddings_list), n_ssr)
        cause_var_embedding_list <- matrix(potential_embeddings_list[embedding_idx,], ncol = ncol(potential_embeddings_list))
      } else {
        cause_var_embedding_list <- matrix(potential_embeddings_list, ncol = ncol(potential_embeddings_list))
      }
    } else {
      cause_var_embedding_list <- NULL
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
    if (is.null(cause_var_embedding_list)) {nrow_res <- 1} else {nrow_res <- nrow(cause_var_embedding_list)}
    multiview_res <- data.frame(matrix(NA, ncol = length(multiview_colnames), nrow = nrow_res))
    colnames(multiview_res) <- multiview_colnames

    ## Main loop of randomized simplex projections
    for (i in 1:nrow_res) {
      ## Choose actual embedding ids (id = 1 is always a target column)
      if (nrow_res > 1) {multiview_idx <- c(1, cause_var_embedding_list[i,])} else {multiview_idx <- 1}
      ## Perform simplex projection
      multiview_res[i,] <- rEDM::block_lnlp(block_multiview[,multiview_idx],
                                            lib = c(1, nrow(block_multiview)),
                                            pred = c(1, nrow(block_multiview)),
                                            method = "simplex",
                                            silent = TRUE,
                                            tp = tp,
                                            target_column = 1)
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
                                    #lib = 1:nrow(block_multiview),
                                    #pred = 1:nrow(block_multiview),
                                    lib = lib,
                                    pred = pred,
                                    target_column = 1,
                                    tp = tp,
                                    dist_w = multiview_dist,
                                    theta = theta,
                                    regularized = regularized,
                                    lambda = lambda,
                                    alpha = alpha,
                                    glmnet_parallel = glmnet_parallel,
                                    save_smap_coefficients = save_smap_coefficients)

    # Prepare recording parameters
    if (is.null(lambda)) lambda <- "NULL"
    if (is.null(alpha)) alpha <- "NULL"
    if (analysis_flow == "mdr") {
      uic_method_record <- "NULL"
    } else {
      uic_method_record <- uic_method
    }
    if (is.null(max_delay)) {
      max_delay_record <- "NULL"
    } else {
      max_delay_record <- max_delay
    }
    # Record parameters
    mdr_res$uic_res <- uic_res
    mdr_res$top_embeddings <- top_multiview_res
    mdr_res$multiview_dist <- multiview_dist
    mdr_res$parms <- data.frame(n_lib = (lib[2]-lib[1]+1),
                                n_pred = (pred[2]-pred[1]+1),
                                E_var = Ex,
                                max_delay = max_delay_record,
                                uic_method = uic_method_record,
                                use_all_lag = use_all_lag,
                                n_ssr = nrow(cause_var_embedding_list),
                                k = nrow(top_multiview_res),
                                evaluate_by = evaluate_by,
                                theta = theta,
                                regularized = regularized,
                                lambda = lambda,
                                alpha = alpha,
                                random_seed = random_seed)

    # Output notes
    if (!silent) message("\nThis function is a beta version and multiple validations are ongoing.")

    # Return results
    return(mdr_res)
  }
}
