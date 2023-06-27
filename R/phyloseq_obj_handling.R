#' @title Bundle taxa names in a phyloseq object
#' @description \code{taxa_name_bundle} bundle taxa names in a phylosq object for a visualization purpose. Minor taxa names will be summarized as "Others". Also, if the taxa name is empty, the function assigns "Undetermined" for the taxa.
#' @importFrom magrittr %>%
#' @param ps_obj Phyloseq object.
#' @param taxa_rank_bundled Character. The name of phylogenetic rank that will be bundled.
#' @param top_taxa_n Numeric. The number of taxa names that will be retained.
#' @param new_taxa_rank Character. The name of the bundled taxa rank. Default is `bundled_tax`.
#' @param taxa_rank_list Character vector. The specified taxa ranks are used to determine whether a tax is "Undetermined" or not.
#' @param add_undet Logical. If TRUE, "Undetermined" is assigned to undetermined taxa. If FALSE, "Others" is assigned.
#' @return Phyloseq object. Bundled taxa names are in the `new_taxa_rank` column.
#' @export
#' @examples
#' # taxa_name_bundle(ps_obj, "phylum", top_taxa_n = 10)
taxa_name_bundle <- function(ps_obj,
                             taxa_rank_bundled,
                             top_taxa_n = 10,
                             new_taxa_rank = "bundled_tax",
                             taxa_rank_list = c("superkingdom",
                                                "kingdom",
                                                "phylum",
                                                "order",
                                                "family",
                                                "genus",
                                                "species"),
                             add_undet = TRUE) {
  tax_df_all <- as.data.frame(phyloseq::tax_table(ps_obj))
  tax_df <- tax_df_all[,taxa_rank_list]
  # Warning message
  message(sprintf("Please ensure that there is no information other than taxa names in the columns between \'%s\' and \'%s\'.", taxa_rank_list[1], rev(taxa_rank_list)[1]))

  # Check rep_tax name
  if (!(new_taxa_rank %in% colnames(tax_df))) {
    if (add_undet) {
      tax_df$tmp_xxxxx <- "Undetermined"
      colnames(tax_df)[ncol(tax_df)] <- new_taxa_rank
      #tax_df <- dplyr::mutate(tax_df, !!new_taxa_rank := "Undetermined")
      others_name <- c("Others", "Undetermined")
    } else {
      tax_df$tmp_xxxxx <- "Others"
      colnames(tax_df)[ncol(tax_df)] <- new_taxa_rank
      #tax_df <- dplyr::mutate(tax_df, !!new_taxa_rank := "Others")
      others_name <- c("Others")
    }
  } else {
    stop("colname \'new_taxa_rank\' is used in taxa_name_bundle(). Please rename the colunm name.")
  }

  # Identify NA or blank elements in taxa_rank_bundled
  rep_tax_cond1 <- tax_df[, taxa_rank_bundled] == "" | is.na(tax_df[, taxa_rank_bundled])
  tax_col1 <- which(colnames(tax_df) == taxa_rank_bundled)
  tax_col2 <- which(colnames(tax_df) == rev(taxa_rank_list)[1])

  # Identify rows of which cells between target and species columns are all NA or blank
  na_cond1 <- tax_df[, (tax_col1+1):tax_col2] == ""
  na_cond2 <- is.na(tax_df[, (tax_col1+1):tax_col2])
  rep_tax_cond2 <- apply(na_cond1 | na_cond2, 1, sum) == (tax_col2 - tax_col1)

  # Add representative taxa information
  tax_df[!rep_tax_cond1, new_taxa_rank] <- as.character(tax_df[!rep_tax_cond1, taxa_rank_bundled])
  tax_df[rep_tax_cond1 & !rep_tax_cond2, new_taxa_rank] <- "Others"

  # Rewrite taxa information
  tax_df_all$tmp_xxxxx <- tax_df[,new_taxa_rank]
  colnames(tax_df_all)[ncol(tax_df_all)] <- new_taxa_rank
  ps_obj2 <- phyloseq::phyloseq(phyloseq::otu_table(ps_obj),
                                phyloseq::sample_data(ps_obj), phyloseq::tax_table(as.matrix(tax_df_all)))

  # Assign "Others" to rare taxa
  taxa_abundance_rank <- stats::aggregate(phyloseq::taxa_sums(ps_obj2),
                                          by = list(phyloseq::tax_table(ps_obj2)[, new_taxa_rank]),
                                          sum)
  taxa_abundance_rank <- taxa_abundance_rank[order(taxa_abundance_rank$x, decreasing = T), ]

  # Check taxa ranks (except "Others" and "Undetermined")
  nd_id <- match(others_name, taxa_abundance_rank$Group.1)
  if (all(is.na(nd_id))) {
    taxa_top <- taxa_abundance_rank[1:top_taxa_n, ]
  } else {
    nd_id <- na.omit(nd_id)
    taxa_top <- taxa_abundance_rank[-nd_id, ][1:top_taxa_n, ]
  }
  taxa_top <- taxa_abundance_rank[-nd_id,][1:top_taxa_n,]
  rare_tax <- is.na(match(phyloseq::tax_table(ps_obj2)[, new_taxa_rank],
                          c(as.character(taxa_top$Group.1), others_name)))
  phyloseq::tax_table(ps_obj2)[rare_tax, new_taxa_rank] <- "Others"

  # Return the phyloseq object
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

