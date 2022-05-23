#' @title Generate twin surrogate time series
#' @description \code{twin_surrogate_cpp} generates twin surrogate time series.
#' @param original Original time series.
#' @param dim Integer. Embedding dimension.
#' @param num.iter Numeric. The number of surrogate to be generated.
#' @param tau Numeric. Time delay in the embedded time series.
#' @param s Numeric. Threshold for the twin in the recurrence plot.
#' @param surrogate.option Character. If `random`, normal twin surrogate will be generated. If `phase_lock`, the periodicity of the time series will be preserved.
#' @param initial.point Character. If `same_season`, twin surrogate time series starts with a point in a same season. If `twins`, twin surrogate time series starts with a twin point.
#' @param distance.method Character. Specify how the distance between points in a state space is calculated.
#' @param point.per.year Numeric. Specify how many observations were made in a unit time (e.g. 24 observations per year).
#' @param s.update Character. If `TRUE`, update `s` if a sufficient number of twins is not found.
#' @param n.twin.threshold Numeric. The minimum number of twins to be accepted.
#' @param output.message Logical. If `TRUE`, show message.
#' @return Twin surrogate time series.
#' @export
#' @examples
#' # twin_surrogate_cpp(time_series, 3, 100)
twin_surrogate_cpp <- function(original,
                               dim,
                               num.iter,
                               tau = 1,
                               s   = 0.875,
                               surrogate.option = c("random", "phase_lock"),
                               initial.point    = c("same_season", "twins"),
                               distance.method  = c("norm", "euclid"),
                               point.per.year   = 24,
                               s.update         = c("on", "off"),
                               n.twin.threshold = 10,
                               output.message = F) {
  surrogate.option = match.arg(surrogate.option)
  initial.point = match.arg(initial.point)
  distance.method = match.arg(distance.method)
  s.update = match.arg(s.update)

  # generate time-lag embedding matrix
  if (dim >  1) {
    original_e <- tseriesChaos::embedd(original, dim, tau)
  }
  if (dim == 1) {
    original_e <- as.matrix(original)
  }
  if (dim <  1) {
    cat("Warning: Embedding dimension should be >= 1")
  }

  # Check arguments
  if ( !(surrogate.option %in% c("random", "phase_lock")) ){
    stop("Error: specify the correct option!")
  }

  if( !(initial.point %in% c("same_season", "twins")) ){
    stop("Error: specify the correct option!")
  }

  surrogate.one.col <- cpp_twinSurrogate(
    original_e,
    dim,
    num.iter,
    s,
    surrogate.option,
    initial.point,
    distance.method,
    point.per.year,
    s.update,
    n.twin.threshold,
    output.message
  )
  return(as.data.frame(surrogate.one.col))
}

