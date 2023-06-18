#' @title Summarize taxa names in a phyloseq object
#' @description \code{taxa_name_summarize} summarizes taxa names in a phylosq object for a visualization purpose. Minor taxa names will be summarized as "Others". Also, if the taxa name is empty, the function assigns "Undetermined" for the taxa.
#' @importFrom magrittr %>%
#' @param ps_obj Phyloseq object.
#' @param taxa_rank Character. The name of phylogenetic rank that will be summarized.
#' @param top_taxa_n Numeric. The number of taxa names that will be retained.
#' @param lowest_rank Character. The name of the lowerst taxa rank. Usually this is `species`.
#' @return Rarefied phyloseq object (`ps_rare `). Summarized taxa names are in the `rep_tax` column.
#' @export
#' @examples
#' # taxa_name_summarize(ps_obj, "phylum", top_taxa_n = 10)
taxa_name_summarize <- function(ps_obj,
                                taxa_rank,
                                top_taxa_n = 10,
                                lowest_rank = "species") {
  tax_df <- as.data.frame(phyloseq::tax_table(ps_obj))

  # Check lowest_rank name
  if (!(lowest_rank %in% colnames(tax_df))) {
    stop(sprintf("column \'%s\' was not found. Please define the colunm name of the lowest rank.", lowest_rank))
  }

  # Check rep_tax name
  if (!("rep_tax" %in% colnames(tax_df))) {
    tax_df$rep_tax <- "Undetermined"
  } else {
    stop("colname \'rep_tax\' is used in taxa_name_summarize(). Please rename the colunm name.")
  }

  # Warning message
  message(sprintf("Please ensure that there is no information other than taxa names in the columns between \'%s\' and \'%s\'.", taxa_rank, lowest_rank))

  rep_tax_cond1 <- tax_df[, taxa_rank] == "" | is.na(tax_df[, taxa_rank])
  tax_col1 <- which(colnames(tax_df) == taxa_rank)
  tax_col2 <- which(colnames(tax_df) == "species")

  # Identify rows of which cells between target and species columns are all NA or blank
  na_cond1 <- tax_df[, (tax_col1+1):tax_col2] == ""
  na_cond2 <- is.na(tax_df[, (tax_col1+1):tax_col2])
  rep_tax_cond2 <- apply(na_cond1 | na_cond2, 1, sum) == (tax_col2 - tax_col1)

  # Add representative taxa information
  tax_df[!rep_tax_cond1, "rep_tax"] <- as.character(tax_df[!rep_tax_cond1, taxa_rank])
  tax_df[rep_tax_cond1 & !rep_tax_cond2, "rep_tax"] <- "Others"
  ps_obj2 <- phyloseq::phyloseq(phyloseq::otu_table(ps_obj),
                                phyloseq::sample_data(ps_obj), phyloseq::tax_table(as.matrix(tax_df)))

  taxa_abundance_rank <- stats::aggregate(phyloseq::taxa_sums(ps_obj2),
                                          by = list(phyloseq::tax_table(ps_obj2)[, "rep_tax"]),
                                          sum)
  taxa_abundance_rank <- taxa_abundance_rank[order(taxa_abundance_rank$x, decreasing = T), ]
  taxa_top <- taxa_abundance_rank[1:top_taxa_n, ]
  low_tax <- is.na(match(phyloseq::tax_table(ps_obj2)[, "rep_tax"],
                         as.character(taxa_top[, 1])))
  phyloseq::tax_table(ps_obj2)[low_tax, "rep_tax"] <- "Others"
  return(ps_obj2)
}



#' @title Assign colors to phyloseq::plot_bar object
#' @description \code{assign_ps_barcol} Assign colors to phyloseq::plot_bar object
#' @param gg_obj ggplot2 object.
#' @param other_name Character. The name of "other" taxa.
#' @param other_col Character. The colors to fill "other" taxa.
#' @param palette Character. Currently only "igv" is available.
#' @param manual_col_code Character vector. The colors to fill "other" taxa.
#' @param legend_title Character. The title of the color legend.
#' @return ggplot2 object.
#' @export
#' @examples
#' # taxa_name_summarize(ps_obj, "phylum", top_taxa_n = 10)
assign_ps_barcol <- function(gg_obj,
                             other_name = c("Others", "Undetermined"),
                             other_col = c("lightgray", "darkgray"),
                             palette = "igv",
                             manual_col_code = NULL,
                             legend_title = NULL) {

  if(length(other_name) != length(other_col)) {
    stop("Lengths of \`other_name\` and \`other_col\` should match!")
  }

  # Collect information
  fill_var_name <- as.character(gg_obj$mapping$fill)[2]
  fill_name <- unique(gg_obj$data[,fill_var_name])
  fill_n <- length(fill_name)
  other_n <- length(other_col)

  # Select colors
  if (!is.null(manual_col_code)) {
    manual_col_code <- manual_col_code[1:fill_n]
    col_code <- manual_col_code
  } else if (palette == "igv") {
    col_code <- ggsci::pal_igv("default")(fill_n)
  } else {
    stop("Specify valid a color palette name")
  }

  # Assign others' colors
  col_code[(fill_n-other_n+1):fill_n] <- other_col

  # Re-order taxa names
  other_id <- match(other_name, fill_name)
  tax_levels <- c(sort(fill_name[-other_id]), fill_name[other_id])
  gg_obj$data[,fill_var_name] <- factor(gg_obj$data[,fill_var_name], levels = tax_levels)

  # Assign colors
  if (is.null(legend_title)) legend_title <- fill_var_name
  gg_obj <- gg_obj + ggplot2::scale_fill_manual(values = col_code, name = legend_title)

  # Return gg_obj
  return(gg_obj)
}

