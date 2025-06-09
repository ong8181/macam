#' @title A custom function for Simplex Projection
#' @description \code{extended_simplex} performs Simplex Projection
#' @param vectors Original time series vector.
#' @param target Target vector.
#' @param lib_indices Library index.
#' @param pred_indices Prediction index.
#' @param num_neighbors The number of nearest neighbors for making predictions.
#' @return A list that contains predictions and statistics
#' @noRd
#' @details
#' \itemize{
#'  \item{Sugihara & May (1990) Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series. Nature 344, 734–741. https://doi.org/10.1038/344734a0}
#' }
#' @keywords internal
extended_simplex <- function(vectors,
                             target,
                             lib_indices,
                             pred_indices,
                             num_neighbors)
{
  # setup output
  pred_vals <- rep.int(NaN, times = length(target))

  # exclude libs that contains NaN in target
  lib_indices <- !apply(cbind(!lib_indices, is.na(target)), 1, any)

  # make predictions
  for (p in which(pred_indices))
  {
    temp_lib <- lib_indices[p]
    lib_indices[p] <- FALSE
    libs <- which(lib_indices)

    # compute distances
    q <- matrix(rep(vectors[p,], length(libs)),
                nrow = length(libs), byrow = T)
    distances <- sqrt(rowSums((vectors[libs,] - q) ^ 2))

    # find nearest neighbors
    neighbors <- order(distances)[1:num_neighbors]
    min_distance <- distances[neighbors[1]]

    # compute weights
    if (!is.na(min_distance) & min_distance == 0) # perfect match
    {
      weights <- rep.int(1e-5, times = num_neighbors)
      weights[distances[neighbors] == 0] <- 1
      # calculate total weights
      total_weight <- sum(weights)
      # make prediction
      pred_vals[p] <- (weights %*% target[libs[neighbors]]) / total_weight
    } else if (!is.na(min_distance) & min_distance > 0) {
      weights <- exp(-distances[neighbors] / min_distance)
      weights[weights < 1e-5] <- 1e-5
      # calculate total weights
      total_weight <- sum(weights)
      # make prediction
      pred_vals[p] <- (weights %*% target[libs[neighbors]]) / total_weight
    } else {
      pred_vals[p] <- NaN
    }

    lib_indices[p] <- temp_lib
  }

  # return output & stats
  return(list(pred = pred_vals,
              stats = compute_stats_SSR(target[pred_indices], pred_vals[pred_indices])))
}
