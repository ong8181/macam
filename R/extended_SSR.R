#' @title Perform regularized S-map using a generalized function
#' @description \code{extended_lnlp} performs the regularized S-map introduced in Censi et al. (2019) Methods in Ecology and Evolution. Multivariate S-map is also supported.
#' @param block_time Dataframe or matrix. Original time series.
#' @param lib Numeric vector. Library indices.
#' @param pred Numeric vector. Prediction indices.
#' @param tp Forecasting time ahead.
#' @param target_column Numeric. Indicates target column
#' @param lib_column Numeric. Indicates library column
#' @param num_neighbors Numeric. The number of nearest neighbors.
#' @param theta Numeric. Weighing function for S-map.
#' @param dist_w Matrix. Distance matrix used to calculate weights for S-map. Implemented for MDR S-map (Chang et al. 2021) Ecology Letters. If `NULL`, then weights are calculated based on Euclidean distance.
#' @param method Character. `s-map` or `simplex`.
#' @param regularized Logical If `TRUE`, regularized S-map will be performed. If `FALSE`, the normal S-map will be performed. Please use `rEDM::s_map` function.
#' @param lambda Numeric. Specify the strength of penalty in the regularization.
#' @param alpha Numeric. `alpha = 0` is the ridge regression, `alpha = 1` is the lasso regression, and `0 < alpha < 1` is an elastic net.
#' @param glmnet_parallel Logical. If TRUE, the computation will be parallel (currently, experimental).
#' @param random_seed Numeric. Random seed.
#' @param save_smap_coefficients Logical. If `TRUE`, S-map coefficients will be saved.
#' @return A list containing:\tabular{ll}{
#'    \code{model_output} \tab  Model predictions \cr
#'    \tab \cr
#'    \code{stats} \tab  Statistics. \cr
#'    \tab \cr
#'    \code{smap_coefficients} \tab  S-map coefficients \cr
#' }
#' @details
#' \itemize{
#'  \item{Cenci, S, Sugihara, G, Saavedra, S. Regularized S-map for inference and forecasting with noisy ecological time series. Methods Ecol Evol. 2019; 10: 650– 660. https://doi.org/10.1111/2041-210X.13150}
#' }
#' @export
#' @examples
#' # extended_lnlp()
extended_lnlp <- function(block_time,
                          lib = c(1, NROW(block_time)),
                          pred = lib,
                          tp = 1,
                          target_column = 1,
                          lib_column = 1:NCOL(block_time),
                          num_neighbors = NCOL(block_time) + 1,
                          theta = 0,
                          dist_w = NULL,
                          method = "s-map", # or "simplex"
                          regularized = FALSE,
                          lambda = NULL,
                          alpha = 0, # default is the ridge regression. If alpha = 1, then do lasso regression
                          glmnet_parallel = FALSE,
                          random_seed = NULL,
                          #no_parallel = glmnet_parallel,
                          save_smap_coefficients = FALSE)
{
  # do multivariate prediction using s-map
  # theta = relative weighting of neighbors based on Euclidean distance in a state space
  #no_parallel = glmnet_parallel

  if(!is.matrix(block_time)) block_time <- as.matrix(block_time)
  n <- NROW(block_time)
  lib <- matrix(lib, ncol = 2)
  pred <- matrix(pred, ncol = 2)

  # setup vectors
  vectors <- matrix(block_time[,lib_column], ncol = length(lib_column))

  # setup lib_indices
  lib_indices <- rep.int(FALSE, times = n)
  for(i in 1:NROW(lib)){
    row_start <- lib[i, 1]
    row_end <- lib[i, 2] - tp
    if(row_end > row_start) lib_indices[row_start:row_end] <- TRUE
  }

  # setup pred_indices
  pred_indices <- rep.int(FALSE, times = n)
  for(i in 1:NROW(pred)){
    row_start <- pred[i, 1]
    row_end <- pred[i, 2] - tp
    if(row_end > row_start) pred_indices[row_start:row_end] <- TRUE
  }

  # setup target
  target <- rep.int(NaN, times = n)
  target[1:(n-tp)] <- block_time[(1+tp):n, target_column]

  if(method == "simplex"){
    simp_out <- extended_simplex(vectors,
                                 target,
                                 lib_indices,
                                 pred_indices,
                                 num_neighbors)
    pred_df <- data.frame(time = 1:nrow(vectors),
                          obs = c(rep(NaN, tp), utils::head(target, n = n - tp)),
                          pred = c(rep(NaN, tp), utils::head(simp_out$pred, n = n - tp)))
    stats <- simp_out$stats

    # Return results
    return(list(model_output = pred_df, stats = stats))

  }else if(method == "s-map" |method == "s_map" |method == "smap"){
    smap_out <- extended_smap(vectors, target,
                              lib_indices, pred_indices, theta,
                              dist_w = dist_w,
                              regularized = regularized,
                              lambda = lambda,
                              alpha = alpha,
                              glmnet_parallel = glmnet_parallel,
                              random_seed = random_seed,
                              #no_parallel = no_parallel,
                              save_smap_coefficients = save_smap_coefficients)
    pred_df <- data.frame(time = 1:nrow(vectors),
                          obs = c(rep(NaN, tp), utils::head(target, n = n - tp)),
                          pred = c(rep(NaN, tp), utils::head(smap_out$pred, n = n - tp)))

    stats <- smap_out$stats

    # Return results
    if(save_smap_coefficients){
      smap_coef_df <- cbind(data.frame(time = 1:nrow(vectors)),
                            matrix(NaN, ncol = dim(smap_out$smap_coefficients)[2],
                                   nrow = dim(smap_out$smap_coefficients)[1]))
      smap_coef_df[,2:(dim(smap_out$smap_coefficients)[2]+1)] <- smap_out$smap_coefficients
      colnames(smap_coef_df)[2:(dim(smap_out$smap_coefficients)[2]+1)] <- colnames(smap_out$smap_coefficients)

      return(list(model_output = pred_df, stats = stats,
                  smap_coefficients = smap_coef_df))
    }else{
      return(list(model_output = pred_df, stats = stats))
    }
  }
}


#' @title Compute statistics
#' @description \code{compute_stats_SSR} computes some statistics.
#' @param obs Observations.
#' @param pred Predictions.
#' @return A data frame containing:\tabular{ll}{
#'    \code{N} \tab  The number of observation \cr
#'    \tab \cr
#'    \code{rho} \tab  Correlation coefficient \cr
#'    \tab \cr
#'    \code{mae} \tab  Mean absolute error \cr
#'    \tab \cr
#'    \code{rmse} \tab  Root mean square error \cr
#' }
#' @noRd
#' @keywords internal
compute_stats_SSR <- function(obs, pred)
{
  # computes performance metrics for how well predictions match observations
  # obs = vector of observations
  # pred = vector of prediction

  N = sum(is.finite(obs) & is.finite(pred))
  rho = stats::cor(obs, pred, use = "pairwise.complete.obs")
  mae = mean(abs(obs-pred), na.rm = TRUE)
  rmse = sqrt(mean((obs-pred)^2, na.rm = TRUE))
  return(data.frame(N = N, rho = rho, mae = mae, rmse = rmse))
}



#' @title Perform regularized S-map
#' @description \code{s_map_rgl} is a wrapper function of `extended_lnlp()` for regularized S-map. For multivariate S-map, please use `extended_lnlp()`.
#' @param ts_data Data.frame, matrix, or vector. One-column univariate time series.
#' @param E Integer. Embedding dimension.
#' @param lib Library indices.
#' @param pred Prediction indices.
#' @param tp Forecasting time ahead.
#' @param num_neighbors Numeric. The number of nearest neighbors.
#' @param theta Numeric. Weighing function for S-map.
#' @param dist_w Matrix. Distance matrix used to calculate weights for S-map. Implemented for MDR S-map (Chang et al. 2021) Ecology Letters. If `NULL`, then weights are calculated based on Euclidean distance.
#' @param regularized Logical If `TRUE`, regularized S-map will be performed. If `FALSE`, the normal S-map will be performed. Please use `rEDM::s_map` function.
#' @param lambda Numeric. Specify the strength of penalty in the regularization.
#' @param alpha Numeric. `alpha = 0` is the ridge regression, `alpha = 1` is the lasso regression, and `0 < alpha < 1` is an elastic net.
#' @param glmnet_parallel Logical. If TRUE, the computation will be parallel (currently, experimental).
#' @param random_seed Numeric. Random seed.
#' @param save_smap_coefficients Logical. If `TRUE`, S-map coefficients will be saved.
#' @return A list containing:\tabular{ll}{
#'    \code{model_output} \tab  Model predictions \cr
#'    \tab \cr
#'    \code{stats} \tab  Statistics. \cr
#'    \tab \cr
#'    \code{smap_coefficients} \tab  S-map coefficients \cr
#' }
#' @export
#' @examples
#' # s_map_rgl()
s_map_rgl <- function(ts_data,
                      E,
                      lib = c(1, length(ts_data)),
                      pred = lib,
                      tp = 1,
                      num_neighbors = length(ts_data) + 1,
                      theta = 0,
                      dist_w = NULL,
                      regularized = FALSE,
                      lambda = NULL,
                      alpha = 0, # default is the ridge regression. If alpha = 1, then do lasso regression
                      glmnet_parallel = FALSE,
                      random_seed = NULL,
                      save_smap_coefficients = FALSE)
{
  # Embed time series
  if (E == 1) {
    ts_embed <- ts_data
  } else if (E > 1) {
    # Convert to  data.frame
    ts_embed <- as.data.frame(ts_data)
    colnames(ts_embed) <- "col1"
    # Add time-delayed time series as a new column
    for(i in 1:(E-1)) {
      ts_embed <- dplyr::mutate(ts_embed, new_col = dplyr::lag(ts_data, n = i))
      colnames(ts_embed)[i+1] <- sprintf("col1_%s", i)
    }
  } else {
    stop("E should be > 0.")
  }

  # Perform S-map
  smap_res <- extended_lnlp(ts_embed,
                            lib = lib,
                            pred = pred,
                            tp = tp,
                            target_column = 1,
                            lib_column = 1:NCOL(ts_embed),
                            num_neighbors = num_neighbors,
                            dist_w = dist_w,
                            theta = theta,
                            regularized = regularized,
                            lambda = lambda,
                            alpha = alpha,
                            glmnet_parallel = glmnet_parallel,
                            random_seed = random_seed,
                            save_smap_coefficients = save_smap_coefficients)

  # Return the result
  return(smap_res)
}

