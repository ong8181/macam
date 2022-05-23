#' @title Power label function for ggplot2
#' @description \code{label_10_to_power} transforms a numeric variable to a scientific notation.
#' @param x Numeric. A numeric variable that will be converted to a scientific notation.
#' @return Scientific notation.
#' @export
#' @examples
#' label_10_to_power(100)
label_10_to_power <- function(x) {
  ifelse(x == 0, "0",
         parse(text = gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))))
  )
}
