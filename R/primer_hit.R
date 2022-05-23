#' @title AllOrients
#' @description \code{AllOrients} expands a primer sequence and show all orientations.
#' @param primer Character. Primer sequences (e.g., ATGC).
#' @export
#' @return An output contains:\tabular{ll}{
#'    \code{Forward} \tab  Forward sequence \cr
#'    \tab \cr
#'    \code{Complement} \tab  Complement sequence \cr
#'    \tab \cr
#'    \code{Reverse} \tab  Reverse sequence \cr
#'    \tab \cr
#'    \code{RevComp} \tab  Reverse-complement sequence \cr
#' }
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
#' @description \code{PrimerHits} search primer sequences in a FASTA/FASTQ file.
#' @param primer Primer sequence
#' @param fn File path to a FASTA/FASTQ file
#' @return The number of primer hits.
#' @export
#' @examples
#' # PrimerHits(primer, fn)
# From https://benjjneb.github.io/dada2/ITS_workflow.html
PrimerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
