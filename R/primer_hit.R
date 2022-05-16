#' @title AllOrients
#' @description \code{AllOrients} xxx
#' @param primer xxx.
#' @export
#' @examples
#' # AllOrients(primer)
# From https://benjjneb.github.io/dada2/ITS_workflow.html
AllOrients <- function(primer) {
  # Create all orientations of the input sequence
  dna <- Biostrings::DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, Biostrings::toString))  # Convert back to character vector
}


#' @title PrimerHits
#' @description \code{PrimerHits} xxx
#' @param primer xxx.
#' @param fn xxx.
#' @export
#' @examples
#' # PrimerHits(primer, fn)
# From https://benjjneb.github.io/dada2/ITS_workflow.html
PrimerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
