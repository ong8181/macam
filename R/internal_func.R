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
