#' Chat.Ind from iNEXT:::Chat.Ind
#' @param x x
#' @param m m
#' @noRd
#' @keywords internal

Chat.Ind <- function (x, m)
{
  x <- x[x > 0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n -
                                                            1)/n * f1^2/2/f2)
  A <- ifelse(f1 > 0, n * f0.hat/(n * f0.hat + f1), 1)
  Sub <- function(m) {
    if (m < n) {
      xx <- x[(n - x) >= m]
      out <- 1 - sum(xx/n * exp(lgamma(n - xx + 1) - lgamma(n -
                                                              xx - m + 1) - lgamma(n) + lgamma(n - m)))
    }
    if (m == n)
      out <- 1 - f1/n * A
    if (m > n)
      out <- 1 - f1/n * A^(m - n + 1)
    out
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
