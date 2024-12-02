#' @title Calculate coverage information of a phyloseq object
#' @description \code{coverage_info2} Calculate coverage information of a phyloseq object
#' @param ps_obj Phyloseq object.
#' @return data.frame or list.
#' @export
#' @examples
#' # coverage_info2(ps_obj)
coverage_info2 <- function (ps_obj) {
  if (phyloseq::taxa_are_rows(ps_obj)) {
    otu_df <- data.frame(t(phyloseq::otu_table(ps_obj)))
  } else {
    otu_df <- data.frame(phyloseq::otu_table(ps_obj))
  }

  # Prepare output object
  cov_df <- data.frame(Assemblage = phyloseq::sample_names(ps_obj),
                       n = phyloseq::sample_sums(ps_obj),
                       S.obs = rowSums(otu_df > 0),
                       f1 = rowSums(otu_df == 1),
                       f2 = rowSums(otu_df == 2),
                       f3 = rowSums(otu_df == 3),
                       f4 = rowSums(otu_df == 4),
                       f5 = rowSums(otu_df == 5),
                       f6 = rowSums(otu_df == 6),
                       f7 = rowSums(otu_df == 7),
                       f8 = rowSums(otu_df == 8),
                       f9 = rowSums(otu_df == 9),
                       f10 = rowSums(otu_df == 10),
                       SC = NA)

  # Add sample coverage
  for (i in 1:nrow(cov_df)) {
    if(cov_df$n[i] > 0) {
      cov_df$SC[i] <- 1 - suppressWarnings(vegan::rareslope(otu_df[i,], cov_df$n[i] - 1))
    }
  }

  return(cov_df)
}

#' @title Perform coverage-based rarefaction for a phyloseq object using vegan
#' @description \code{rarefy_even_coverage2} performs coverage-based rarefaction for a phyloseq object using vegan package
#' @importFrom magrittr %>%
#' @param ps_obj Phyloseq object.
#' @param coverage Numeric. Coverage specified by a user (default = 0.97 = 97%).
#' @param remove_not_rarefied Logical. If `TRUE`, samples of which coverage is lower than the specified coverage will be removed.
#' @param include_rarecurve_results Logical. If TURE, `rarecurve` results will also be included for visualization. The first object is a rarefied phyloseq object, and the second object is `rarecurve` results. If `FALSE`, it returns a rarefied phyloseq object only.
#' @param knots Numeric. Specify `knots` of `rareslope` function.
#' @param n_rarefy_iter Numeric. The number of iterations of rarefactions (default = 1).
#' @param rarefy_average_method  Character. If `n_rarefy_iter` >= 2, this argument determines how the multiple rarefactions are summarized. r`arefy_average_method = "round"` uses `round()`. `rarefy_average_method = "floor"` uses `floor()`. `rarefy_average_method = "ceiling"` uses `ceiling()`.
#' @param rareplot_step_size  Numeric. Step size for `rarecurve` function. Only affect the resolution of a rarefaction plot.
#' @param ran_seed Numeric. Random seed.
#' @return Rarefied phyloseq object (`ps_rare`). If `include_rarecurve_results = TRUE`, `rarecurve` results are stored in the second element of the list.
#' @export
#' @examples
#' # rarefy_even_coverage(ps_obj, coverage = 0.97)
rarefy_even_coverage2 <- function(ps_obj,
                                  coverage = 0.97,
                                  remove_not_rarefied = FALSE,
                                  include_rarecurve_results = FALSE,
                                  knots = 1000,
                                  n_rarefy_iter = 10,
                                  rareplot_step_size = 100,
                                  rarefy_average_method = "round",
                                  ran_seed = 1234) {
  if (!(remove_not_rarefied %in% c(TRUE, FALSE))) {
    stop("Invalid 'remove_not_rarefied'. Should be TRUE or FALSE.")
  }
  if (!(include_rarecurve_results %in% c(TRUE, FALSE))) {
    stop("Invalid 'include_rarecurve_results'. Should be TRUE or FALSE.")
  }

  # Set random seed
  set.seed(ran_seed)

  # Check phyloseq objects
  if (!phyloseq::taxa_are_rows(ps_obj)) {
    com_mat <- ps_obj %>% phyloseq::otu_table() %>% data.frame
  } else {
    com_mat <- ps_obj %>% phyloseq::otu_table() %>% t %>% data.frame
  }

  sam_tbl <- phyloseq::sample_data(ps_obj) %>% as.data.frame()
  all_names <- phyloseq::sample_names(ps_obj)

  if (any(!phyloseq::sample_sums(ps_obj) > 0)) {
    # Store the original object with "zero" samples
    com_mat0 <- com_mat; sam_tbl0 <- sam_tbl; all_names0 <- all_names
    # Extract non-zero samples
    non_zero_sample <- rownames(com_mat[rowSums(com_mat) > 0,])
    com_mat <- com_mat[non_zero_sample,]
    sam_tbl <- sam_tbl[non_zero_sample,]
    all_names <- all_names[rowSums(com_mat0) > 0]
    message("Some samples contain no reads. These samples will NOT be rarefied, but the sample data will be retained in the output if `remove_not_rarefied = FALSE`.\n")
  }

  # Calculate slope
  rare_knots <- knots # knots for slope check
  rareslopelist<-list()
  for(i in 1:nrow(com_mat)){
    rare_p <- round(seq(1, sum(com_mat[i,]) - 1, length.out = rare_knots))
    rareslopelist[[i]] <- data.frame(n_reads = rare_p,
                                     slope = suppressWarnings(vegan::rareslope(com_mat[i,], rare_p)))
    rm(rare_p)
  }

  # Specify coverage
  cvr <- 1 - coverage # Rename the object
  # Identify sequence read number that achieve the coverage
  cvrfun <- function(x) min(which(x$slope < cvr))
  cvrrare <- suppressWarnings(unlist(lapply(rareslopelist, cvrfun)))

  # Summarize sample information
  cvr_df <- data.frame(sample = rownames(com_mat),
                       n_reads_id = cvrrare,
                       rarefy = is.finite(cvrrare),
                       n_reads = NA,
                       coverage = NA)
  rarefy_id <- is.finite(cvrrare)

  # Add sample reads
  for (i in 1:length(cvrrare)) {
    if (is.finite(cvrrare[i])) {
      cvr_df$n_reads[i] <- rareslopelist[[i]]$n_reads[cvrrare[i]]
      cvr_df$coverage[i] <- 1 - suppressWarnings(vegan::rareslope(com_mat[i,], cvr_df$n_reads[i]))
    } else {
      cvr_df$n_reads[i] <- utils::tail(rareslopelist[[i]]$n_reads, n = 1)
      cvr_df$coverage[i] <- 1 - suppressWarnings(vegan::rareslope(com_mat[i,], cvr_df$n_reads[i]))
    }
  }

  # Get rarefied OTU table
  rarefied_df <- data.frame()
  for(i in 1:nrow(cvr_df)) {
    rarefied_df <- rbind(rarefied_df, suppressWarnings(vegan::rrarefy(com_mat[i,], cvr_df$n_reads[i])))
  }

  if (n_rarefy_iter > 1) {
    for (k in 1:n_rarefy_iter) {
      rarefied_df_tmp <- data.frame()
      for(i in 1:nrow(cvr_df)) {
        rarefied_df_tmp <- rbind(rarefied_df_tmp, suppressWarnings(vegan::rrarefy(com_mat[i,], cvr_df$n_reads[i])))
      }
      rarefied_df <- rarefied_df + rarefied_df_tmp
    }
    # Calculate average values
    if (rarefy_average_method == "round") {
      rarefied_df <- round(rarefied_df/n_rarefy_iter)
    } else if (rarefy_average_method == "ceiling") {
      rarefied_df <- ceiling(rarefied_df/n_rarefy_iter)
    } else if (rarefy_average_method == "floor") {
      rarefied_df <- floor(rarefied_df/n_rarefy_iter)
    }
  }

  # Add sample information
  sam_tbl$original_reads <- rowSums(com_mat)
  sam_tbl$original_n_taxa <- rowSums(com_mat > 0)
  for (i in 1:nrow(cvr_df)) sam_tbl$original_coverage[i] <- 1 - suppressWarnings(vegan::rareslope(com_mat[i,], sum(com_mat[i,]) - 1))
  sam_tbl$rarefied <- cvr_df$rarefy
  sam_tbl$rarefied_reads <- cvr_df$n_reads
  sam_tbl$rarefied_n_taxa <- rowSums(rarefied_df > 0)
  sam_tbl$rarefied_coverage <- cvr_df$coverage

  if (any(!phyloseq::sample_sums(ps_obj) > 0)) {
    # Curate data of the original OTU table
    sam_tbl0$original_reads <- rowSums(com_mat0)
    sam_tbl0$original_n_taxa <- rowSums(com_mat0 > 0)
    sam_tbl0$original_coverage <- NA
    # Replace elements of non-zero samples in the original objects
    com_mat0[rownames(rarefied_df),] <- rarefied_df
    sam_tbl0$rarefied <- FALSE
    sam_tbl0$rarefied_reads <- 0
    sam_tbl0$rarefied_n_taxa <- rowSums(com_mat0 > 0)
    sam_tbl0$rarefied_coverage <- NA
    sam_tbl0[rownames(rarefied_df),]$original_coverage <- sam_tbl$original_coverage
    sam_tbl0[rownames(rarefied_df),]$rarefied <- cvr_df$rarefy
    sam_tbl0[rownames(rarefied_df),]$rarefied_reads <- cvr_df$n_reads
    sam_tbl0[rownames(rarefied_df),]$rarefied_coverage <- cvr_df$coverage

    # Replace data in the phyloseq object
    sam_tbl <- sam_tbl0
    rarefied_df <- com_mat0
  }

  ps_rare <- ps_obj
  phyloseq::sample_data(ps_rare) <- phyloseq::sample_data(sam_tbl)

  # Compile output
  if (phyloseq::taxa_are_rows(ps_obj)) {
    phyloseq::otu_table(ps_rare) <- t(rarefied_df) %>% phyloseq::otu_table(taxa_are_rows = TRUE)
  } else {
    phyloseq::otu_table(ps_rare) <- rarefied_df %>% phyloseq::otu_table(taxa_are_rows = FALSE)
  }

  if (any(!rarefy_id)) {
    message1 <- sprintf("Following samples were NOT rarefied because of the low coverage: %s\n",
                        all_names[!rarefy_id] %>% paste(collapse = " "))
  } else {
    message1 <- "All samples were successfully rarefied!\n"
  }
  message(message1)

  if (remove_not_rarefied) {
    ps_rare <- phyloseq::prune_samples(all_names[rarefy_id], ps_rare)
    message("Not-rarefied samples were removed from output as you specified.")
  } else {
    message2 <- "Rarefied/not-rarefied samples were kept in the phyloseq object."
    message3 <- "Sequence reads of the not-rarefied samples were not changed."
    message4 <- "Please check 'rarefied' column of sample_data()."
    message(message2)
    message(message3)
    message(message4)
  }

  if (include_rarecurve_results) {
    if (remove_not_rarefied) {
      com_mat2 <- com_mat[all_names[rarefy_id],]
      sam_tbl2 <- sam_tbl[all_names[rarefy_id],]
    } else {
      com_mat2 <- com_mat
      sam_tbl2 <- sam_tbl[rownames(com_mat),]
    }

    # Collect information
    rcurve_df <- suppressWarnings(vegan::rarecurve(com_mat2, step = rareplot_step_size, tidy = TRUE))
    rpoint_df <- data.frame(sample = rownames(sam_tbl2),
                            rarefied_slope = 1 - sam_tbl2$rarefied_coverage,
                            rarefied_reads = sam_tbl2$rarefied_reads,
                            rarefied_n_taxa = sam_tbl2$rarefied_n_taxa,
                            predicted_n_taxa = NA)
    for(i in 1:nrow(cvr_df)){
      rpoint_df$predicted_n_taxa[i] <- suppressWarnings(vegan::rarefy(com_mat2[i,], cvr_df$n_reads[i]))
    }
    # Compile phyloseq object
    ps_rare <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_rare) > 0, ps_rare)
    vegan_df <- list(rcurve_df, rpoint_df)
    names(vegan_df) <- c("rcurve_df", "rpoint_df")
    return(list(ps_rare, vegan_df))
  } else {
    # Compile phyloseq object
    ps_rare <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_rare) > 0, ps_rare)
    return(ps_rare)
  }
}


#' @title Visualize rarefaction curve
#' @description \code{plot_rarefy} visualizes a coverage-based rarefaction.
#' @importFrom magrittr %>%
#' @param ps_obj Phyloseq object.
#' @param plot_rarefied_point Logical. Specify whether rarefied reads are plotted.
#' @param plot_slope Logical. Specify whether tangent lines are plotted.
#' @param linetype Numeric or character. Specify the line type of the tangent lines.
#' @return ggplot object (`g_rare`). User may further edit the figure using `ggplot2`.
#' @export
#' @examples
#' # plot_rarefy(ps_obj)
plot_rarefy2 <- function (ps_obj,
                          plot_rarefied_point = TRUE,
                          plot_slope = FALSE,
                          linetype = 1) {
  # Retrieve data for ggplot
  rcurve_df <- ps_obj[[2]]$rcurve_df
  rpoint_df <- ps_obj[[2]]$rpoint_df
  slope_vals <- rpoint_df$rarefied_slope
  intercept_vals <- rpoint_df$predicted_n_taxa - rpoint_df$rarefied_reads * rpoint_df$rarefied_slope

  # Visualize
  g_rare <- rcurve_df %>%
    ggplot2::ggplot(ggplot2::aes(x = Sample, y = Species, color = Site)) + ggplot2::geom_line()

  # Add rarefied points
  if (plot_rarefied_point) {
    g_rare <- g_rare +
      ggplot2::geom_point(data = rpoint_df, ggplot2::aes(x = rarefied_reads, y = predicted_n_taxa),
                          color = "gray30", alpha = 0.6)
  }

  # Add slopes
  if (plot_slope) {
    g_rare <- g_rare +
      ggplot2::geom_abline(slope = slope_vals,
                           intercept = intercept_vals, color = "gray60",
                           linetype = linetype, linewidth = 0.3, alpha = 0.5)
  }

  return(g_rare)
}
