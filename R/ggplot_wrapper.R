#' @title Power label function for ggplot2
#' @description \code{label_to_power10} transforms a numeric variable to a scientific notation.
#' @param x Numeric. A numeric variable that will be converted to a sceitific notation.
#' @return Scientific notation.
#' @export
#' @examples
#' label_to_power10(100)
label_to_power10 <- function(x) {
  ifelse(x == 0, "0",
         parse(text = gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))))
  )
}
