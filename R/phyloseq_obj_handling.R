#' @title Summarize taxa names in a phyloseq object
#' @description \code{taxa_name_summarize} summarizes taxa names in a phylosq object for a visualization purpose. Minor taxa names will be summarized as "Others". Also, if the taxa name is empty, the function assigns "Undetermined" for the taxa.
#' @importFrom magrittr %>%
#' @param ps_obj Phyloseq object.
#' @param taxa_rank Character. The name of phylogenetic rank that will be summarized.
#' @param top_taxa_n Numeric. The number of taxa names that will be ratained.
#' @return Rarefied phyloseq object (`ps_rare `). Summarized taxa names are in the `rep_tax` column.
#' @export
#' @examples
#' # taxa_name_summarize(ps_obj, "phylum", top_taxa_n = 10)
taxa_name_summarize <- function(ps_obj,
                                taxa_rank,
                                top_taxa_n = 10){
  tax_df <- as.data.frame(phyloseq::tax_table(ps_obj))
  if(is.null(tax_df$rep_tax)) tax_df$rep_tax <- "Undetermined"

  # Search Others and Undetermined taxa
  rep_tax_cond1 <- tax_df[,taxa_rank] == "" & !is.na(tax_df[,taxa_rank])
  tax_col1 <- which(colnames(tax_df) == taxa_rank)
  tax_col2 <- which(colnames(tax_df) == "species")
  rep_tax_cond2 <- apply(tax_df[,(tax_col1+1):tax_col2] == "", 1, sum) == (tax_col2 - tax_col1)

  # Replace taxa names
  tax_df[!rep_tax_cond1, "rep_tax"] <- as.character(tax_df[!rep_tax_cond1, taxa_rank])
  tax_df[rep_tax_cond1 & !rep_tax_cond2, "rep_tax"] <- "Others"

  # Re-import phyloseq object with revised tax_table
  ps_obj2 <- phyloseq::phyloseq(phyloseq::otu_table(ps_obj),
                                   phyloseq::sample_data(ps_obj),
                                   phyloseq::tax_table(as.matrix(tax_df)))

  # Repalce low abundance taxa name with Others
  taxa_abundance_rank <- stats::aggregate(phyloseq::taxa_sums(ps_obj2),
                                          by = list(phyloseq::tax_table(ps_obj2)[,"rep_tax"]), sum)
  taxa_abundance_rank <- taxa_abundance_rank[order(taxa_abundance_rank$x, decreasing = T),]
  taxa_top <- taxa_abundance_rank[1:top_taxa_n,]

  low_tax <- is.na(match(phyloseq::tax_table(ps_obj2)[,"rep_tax"], as.character(taxa_top[,1])))
  phyloseq::tax_table(ps_obj2)[low_tax,"rep_tax"] <- "Others"

  return(ps_obj2)
}
