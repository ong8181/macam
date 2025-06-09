#' Chat.Ind from iNEXT:::Chat.Ind
#' @param x Numeric vector. Community data of a single site (e.g., `x = c(20,10,4,2,2,3,13,1,0)`).
#' @param m Numeric. Sample size extracted from `x`.
#' @noRd
#' @keywords internal
Chat.Ind <- function (x, m)
{
  x <- x[x > 0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n-1)/n*f1*(f1-1)/2, (n-1)/n*f1^2/2/f2)
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m) {
    if (m < n) {
      xx <- x[(n-x) >= m]
      out <- 1 - sum(xx/n*exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if (m == n) {out <- 1-f1/n*A}
    if (m > n) {out <- 1-f1/n*A^(m-n+1)}
    return(out)
  }
  sapply(m, Sub)
}



#' rarefaction_subsample from phyloseq:::rarefaction_subsample
#' @param x x
#' @param sample.size sample.size
#' @param replace replace
#' @noRd
#' @keywords internal
rarefaction_subsample <- function (x, sample.size, replace = FALSE)
{
  rarvec <- numeric(length(x))
  if (sum(x) <= 0) {
    return(rarvec)
  }
  if (replace) {
    suppressWarnings(subsample <- sample(1:length(x), sample.size,
                                         replace = TRUE, prob = x))
  }
  else {
    obsvec <- apply(data.frame(OTUi = 1:length(x), times = x),
                    1, function(x) {
                      rep_len(x["OTUi"], x["times"])
                    })
    obsvec <- unlist(obsvec, use.names = FALSE)
    suppressWarnings(subsample <- sample(obsvec, sample.size,
                                         replace = FALSE))
  }
  sstab <- table(subsample)
  rarvec[methods::as(names(sstab), "integer")] <- sstab
  return(rarvec)
}


#' singleton_estimator is a function of obtaining the estimator of singleton based the method of Chiu & Chao (2016)
#' @param x a vector of species sample frequencies
#' @param datatype Character. Specify the data type. `abundance` or `incidence` (this should be normally `abundance` for a phyloseq object).
#' @return singleton estimate and a numerical vector of corrected data.
#' @details
#' \itemize{
#'  \item{Chiu & Chao (2016) PeerJ https://doi.org/10.7717/peerj.1634}
#' }
singleton_estimator <- function(x, datatype = "abundance")
{
  # Check data format
  if (is.matrix(x) == T || is.data.frame(x) == T){
    if (ncol(x) != 1 & nrow(x) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(x) == 1) {x <- x[, 1]} else {x <- x[1, ]}
  }
  x <- as.numeric(x)

  # Correct singlton according to Chiu & Chao (2016)
  if(datatype == "abundance"){
    f2 <- sum(x == 2) # Doubleton
    f3 <- sum(x == 3)
    f4 <- sum(x == 4)
    f1 = ifelse(f3*f4>0, 4*f2^2/(3*f3)-f2*f3/(2*f4), 4*f2^2/(3*(f3+1))-f2*f3/(2*(f4+1)))
    I = which(x == 1)
    if(length(I) > 0) x <- x[-I]
    if(f1 > 0) x <- c(x, rep(1, round(f1)))
    x = x[x > 0]
    x <- as.integer(x)
  }else{
    N = x[1]; x = x[-1];
    f2 <- sum(x == 2) # Doubleton
    f3 <- sum(x == 3)
    f4 <- sum(x == 4)
    f1 = ifelse(f3*f4>0, 4*f2^2/(3*f3)-f2*f3/(2*f4), 4*f2^2/(3*(f3+1))-f2*f3/(2*(f4+1)))
    I = which(x == 1)
    if(length(I) > 0) x <- x[-I]
    if(f1 > 0) x <- c(x, rep(1, round(f1)))
    x = x[x > 0]
    x <- as.integer(x)
    x = c(N,x)
  }
  return(list(singleton_est = f1, corrected_data = x))
}


#' @title `make_block` function from rEDM v0.7.5
#' @description \code{make_block} generates a lagged block with the appropriate max_lag and tau, while respecting lib (by inserting NANs, when trying to lag past lib regions)
#' @param block a data.frame or matrix where each column is a time series
#' @param t Numeric. The time index for the block.
#' @param max_lag The total number of lags to include for each variable.
#' @param tau The lag to use for time delay embedding.
#' @param lib A 2-column matrix (or 2-element vector) where each row specifies the first and last *rows* of the time series to use for attractor reconstruction.
#' @param restrict_to_lib Whether to restrict the final lagged block to just the rows specified in lib (if lib exists).
#' @return A data.frame with the lagged columns and a time column. If the original block had columns X, Y, Z and max_lag = 3, then the returned data.frame will have columns TIME, X, X_1, X_2, Y, Y_1, Y_2, Z, Z_1, Z_2.
#' @details
#' \itemize{
#'  \item{Ye et al. (2018) Applications of Empirical Dynamic Modeling from Time Series. https://github.com/ha0ye/rEDM}
#' }
#' @keywords internal
make_block <- function (block, t = NULL, max_lag = 3, tau = 1, lib = NULL,
                        restrict_to_lib = TRUE) {
  if (is.vector(block)) block <- matrix(block, ncol = 1)
  num_vars <- NCOL(block)
  num_rows <- NROW(block)

  if (!is.null(lib)) {
    if (is.vector(lib)) lib <- matrix(lib, ncol = 2, byrow = TRUE)
  }
  output <- matrix(NA, nrow = num_rows, ncol = 1 + num_vars * max_lag)
  col_names <- character(1 + num_vars * max_lag)
  if (is.null(t)) output[, 1] <- 1:num_rows else output[, 1] <- t
  col_names[1] <- "time"
  col_index <- 2
  if (is.null(colnames(block))) colnames(block) <- paste0("col", seq_len(num_vars))

  for (j in 1:num_vars) {
    ts <- block[, j]
    if (is.list(ts)) {
      ts <- unlist(ts)
    }
    output[, col_index] <- ts
    col_names[col_index] <- colnames(block)[j]
    col_index <- col_index + 1
    if (max_lag > 1) {
      for (i in 1:(max_lag - 1)) {
        ts <- c(rep_len(NA, tau), ts[1:(num_rows - tau)])
        if (!is.null(lib)) {
          for (k in seq_len(NROW(lib))) {
            ts[lib[k, 1] - 1 + (1:tau)] <- NA
          }
        }
        output[, col_index] <- ts
        col_names[col_index] <- paste0(colnames(block)[j], "_", i * tau)
        col_index <- col_index + 1
      }
    }
  }
  if (!is.null(lib) && restrict_to_lib) {
    row_idx <- sort(unique(do.call(c, mapply(seq, lib[, 1], lib[, 2], SIMPLIFY = FALSE))))
    output <- output[row_idx, ]
  }
  output <- data.frame(output)
  names(output) <- col_names
  return(output)
}
