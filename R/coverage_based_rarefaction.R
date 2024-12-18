#' @title Calculate coverage information of a phyloseq object
#' @description \code{coverage_info} a wrapper function of `iNEXT::DataInfo()` for a phyloseq object
#' @param ps_obj Phyloseq object.
#' @param datatype Character. Specify the data type. `abundance` or `incidence` (this should be normally `abundance` for a phyloseq object).
#' @param estimate_div logical. If `TRUE`, return the diversity estimate at the minimum coverage using iNEXT::estimateD.
#' @return data.frame or list.
#' @details
#' \itemize{
#'  \item{Hsieh et al. (2016) Methods in Ecology and Evolution https://doi.org/10.1111/2041-210X.12613}
#' }
#' @export
#' @examples
#' # coverage_info(ps_obj)
## the code is based on the iNEXT package
## https://github.com/JohnsonHsieh/iNEXT/blob/master/R/EstIndex.R
coverage_info <- function(ps_obj, datatype = "abundance", estimate_div = FALSE){
  # Exporting otu_table
  if(phyloseq::taxa_are_rows(ps_obj)) {
    # Taxa are in rows
    otu_df <- data.frame(phyloseq::otu_table(ps_obj))
  } else {
    # Taxa are in columns
    otu_df <- data.frame(t(phyloseq::otu_table(ps_obj)))
  }
  # Converting dataframe to list
  otu_list <- otu_df %>% purrr::array_tree(2)
  # iNEXT::DataInfo()
  cov_info <- iNEXT::DataInfo(otu_list, datatype = datatype)
  # Correct singltons
  corrected_otu_list <- purrr::map(otu_list, function(x) singleton_estimator(x)$corrected_data)
  cov_info2 <- iNEXT::DataInfo(corrected_otu_list, datatype = datatype)
  cov_info$f1_cor <- cov_info2$f1
  cov_info$SC_cor <- cov_info2$SC

  if(estimate_div) {
    # Estimate diversity
    estimate_div_df <- iNEXT::estimateD(otu_list, datatype = datatype, base = "coverage",
                                        q = 0, conf = 0.95, level = min(cov_info$SC))
    # Return statistics
    cov_info <- list(cov_info, estimate_div_df)
    names(cov_list) <- c("cov_info", "estimate_div")
    return(cov_list)
  } else {
    # Return statistics
    return(cov_info)
  }
}


#' @title Perform coverage-based rarefaction for a phyloseq object
#' @description \code{rarefy_even_coverage} performs coverage-based rarefaction for a phyloseq object using iNEXT package
#' @importFrom magrittr %>%
#' @param ps_obj Phyloseq object.
#' @param coverage Numeric. Coverage specified by a user (default = 0.97 = 97%).
#' @param remove_not_rarefied Logical. If `TRUE`, samples of which coverage is lower than the specified coverage will be removed.
#' @param correct_singletons Logical. If TURE, singleton counts will be corrected according to the method of Chiu & Chao (2016). Corrected coverage and reads will be used for the rarefaction of the original data.
#' @param include_iNEXT_results Logical. If TURE, iNEXT results will be included as a element of the output. The first object is a rarefied phyloseq object, and the second object is an iNEXT result. Also, if `TRUE`, computation time will increase. If `FALSE`, it returns a rarefied phyloseq object only.
#' @param se Logical. Specify whether SE of the rarefaction curve is calculated.
#' @param nboot Numeric. Specify `nboot` of iNEXT function (valid only if `include_iNEXT_results = TRUE`).
#' @param knots Numeric. Specify `knots` of iNEXT function (valid only if `include_iNEXT_results = TRUE`).
#' @param n_rarefy_iter Numeric. The number of iterations of rarefactions (default = 1).
#' @param rarefy_average_method  Character. If `n_rarefy_iter` >= 2, this argument determines how the multiple rarefactions are summarized. r`arefy_average_method = "round"` uses `round()`. `rarefy_average_method = "floor"` uses `floor()`. `rarefy_average_method = "ceiling"` uses `ceiling()`.
#' @param sample_method Character. Specify which function is used for rarefaction. `sample_method = "vegan"` uses `vegan::rrarefy()`, while `sample_method = "phyloseq"` uses `phyloseq:::rarefaction_subsample()`.
#' @param ran_seed Numeric. Random seed.
#' @return Rarefied phyloseq object (`ps_rare`). If `include_iNEXT_results = TRUE`, `iNEXT` results are stored in the second element of the list.
#' @export
#' @details
#' \itemize{
#'  \item{Hsieh et al. (2016) Methods in Ecology and Evolution https://doi.org/10.1111/2041-210X.12613}
#' }
#' @examples
#' # rarefy_even_coverage(ps_obj, coverage = 0.97)
rarefy_even_coverage <-  function(ps_obj,
                                  coverage = 0.97,
                                  remove_not_rarefied = FALSE,
                                  include_iNEXT_results = FALSE,
                                  correct_singletons = FALSE,
                                  se = FALSE,
                                  nboot = 40,       # Only valid if include_rarefaction_curve = TRUE
                                  knots = 50,       # Only valid if include_rarefaction_curve = TRUE
                                  n_rarefy_iter = 1,
                                  rarefy_average_method = "round",
                                  sample_method = "phyloseq",
                                  ran_seed = 1234
){
  # Check arguments
  if (!(remove_not_rarefied %in% c(TRUE, FALSE))) {
    stop("Invalid \'remove_not_rarefied\'. Should be TRUE or FALSE.")
  }
  if (!(include_iNEXT_results %in% c(TRUE, FALSE))) {
    stop("Invalid \'include_iNEXT_results\'. Should be TRUE or FALSE.")
  }
  if (!(sample_method %in% c("rrarefy", "phyloseq"))) {
    stop("Invalid \'sample_method\'. Should be \'rrarefy\' or \'phyloseq\'.")
  }
  if (!(rarefy_average_method %in% c("round", "ceiling", "floor"))) {
    stop("Invalid \'rarefy_average_method\'. Should be \'round\', \'ceiling\', or \'floor\'.")
  }

  # Set random seed
  set.seed(ran_seed)

  # Convert ps_obj to OTU table
  ## Check taxa_are_rows
  if (phyloseq::taxa_are_rows(ps_obj)) {
    com_mat <- ps_obj %>% phyloseq::otu_table() %>% data.frame
  } else {
    com_mat <- ps_obj %>% phyloseq::otu_table() %>% t %>% data.frame
  }
  all_names <- colnames(com_mat)

  # Get community list
  com_list <- com_mat %>% purrr::array_tree(2)

  # Correct singleton
  if (correct_singletons) {
    com_list <- purrr::map(com_list, function(x) singleton_estimator(x)$corrected_data)
  }

  # Check maximum coverage
  inext_max_sc <- com_list %>% purrr::map(function(x) Chat.Ind(x, sum(x))) %>% unlist
  # Check whether (specified coverage) < (max coverage)
  rarefy_id <- (coverage < inext_max_sc)
  if (all(!rarefy_id)) {
    stop("Depths of all the samples were not sufficient for the rarefaction! Try a decreased coverage.")
  }

  # Do iNEXT
  if (include_iNEXT_results) {
    inext_out <- com_list %>%
      purrr::map(function(x) iNEXT::iNEXT(x, q = 0,
                                          datatype="abundance",
                                          endpoint = sum(x),
                                          se = se,
                                          #conf = 0.95,
                                          nboot = nboot,
                                          knots = knots))
  }
  ## Calculate sample size to archive the specified coverage
  inext_reads <- com_list[rarefy_id] %>%
    purrr::map(function(x) coverage_to_samplesize(x, coverage)) %>% unlist
  inext_coverage <- purrr::map2(com_list[rarefy_id], inext_reads, function(x,y) Chat.Ind(x,y)) %>% unlist

  # Get rarefied counts
  # ! Use original com_mat (NOT com_list).
  # ! If correct_singletons = TRUE,
  # ! then we use "corrected" inext_reads and inext_coverage for the original com_mat
  # ! Estimating which singleton OTU is "true" singlton is currently NOT possible.
  rrlist <- list(x = com_mat[,rarefy_id] %>% t %>% purrr::array_tree(1),
                 y = inext_reads %>% purrr::array_tree())
  ## This increases computation time
  if (sample_method == "rrarefy") {
    for (j in 1:n_rarefy_iter) {
      rarefied_count_list_tmp <- rrlist %>% purrr::pmap(function(x,y) vegan::rrarefy(x, y))
      if (j == 1) {
        rarefied_count_list <- rarefied_count_list_tmp # temporal object
      } else {
        for (k in 1:length(rrlist[[1]])) {
          rarefied_count_list[[k]] <- rarefied_count_list[[k]] %>%
            rbind(rarefied_count_list_tmp[[k]]) %>% colSums
        }
      }
    }
  } else if (sample_method == "phyloseq") {
    for (j in 1:n_rarefy_iter) {
      rarefied_count_list_tmp <- rrlist %>%
        purrr::pmap(function(x,y) rarefaction_subsample(x, y, replace = FALSE))
      if (j == 1) {
        rarefied_count_list <- rarefied_count_list_tmp # temporal object
      } else {
        for (k in 1:length(rrlist[[1]])) {
          rarefied_count_list[[k]] <- rarefied_count_list[[k]] %>%
            rbind(rarefied_count_list_tmp[[k]]) %>% colSums
        }
      }
    }
  }

  ## Take ceiling of the average
  if (rarefy_average_method == "round") {
    rarefied_count_list <- rarefied_count_list %>% purrr::map(function(x) round(x/n_rarefy_iter))
  } else if (rarefy_average_method == "ceiling") {
    rarefied_count_list <- rarefied_count_list %>% purrr::map(function(x) ceiling(x/n_rarefy_iter))
  } else if (rarefy_average_method == "floor") {
    rarefied_count_list <- rarefied_count_list %>% purrr::map(function(x) floor(x/n_rarefy_iter))
  }
  ## Other options

  # Combined as data.frame
  rarefied_count <- as.data.frame(do.call(rbind, rarefied_count_list))
  rownames(rarefied_count) <- names(inext_reads)
  rarefied_sp <- rowSums(rarefied_count > 0)

  # Replace original data
  com_mat2 <- com_mat
  rarefied_count_t <- t(as.matrix(rarefied_count))
  com_mat2[,colnames(rarefied_count_t)] <- rarefied_count_t
  # Import to the original phyloseq object
  ps_rare <- ps_obj
  if (phyloseq::taxa_are_rows(ps_obj)) {
    phyloseq::otu_table(ps_rare) <- com_mat2 %>% phyloseq::otu_table(taxa_are_rows = TRUE)
  } else {
    phyloseq::otu_table(ps_rare) <- t(com_mat2) %>% phyloseq::otu_table(taxa_are_rows = FALSE)
  }
  # Add rarefied/not-rarefied information
  phyloseq::sample_data(ps_rare) <- phyloseq::sample_data(ps_rare) %>%
    data.frame %>%
    dplyr::mutate(rarefied = as.logical(rarefy_id),
                  original_reads = phyloseq::sample_sums(ps_obj),
                  original_n_taxa = rowSums(phyloseq::otu_table(ps_obj) > 0),
                  original_coverage = inext_max_sc,
                  rarefied_reads = phyloseq::sample_sums(ps_rare),
                  rarefied_n_taxa = rowSums(phyloseq::otu_table(ps_rare) > 0),
                  rarefied_pred_n_taxa = rrlist %>% purrr::pmap(function(x,y) suppressWarnings(vegan::rarefy(x, y))) %>% unlist,
                  rarefied_coverage = NA,
                  correct_singletons = correct_singletons)
  phyloseq::sample_data(ps_rare)[rarefy_id,"rarefied_coverage"] <- inext_coverage

  # Output message
  if (any(!rarefy_id)) {
    message1 <- sprintf("Following samples were NOT rarefied: %s\n",
                        all_names[!rarefy_id] %>%
                          paste(collapse=" "))
  } else {
    message1 <- "All samples were successfully rarefied!\n"
  }
  # Show message
  message(message1)

  # Add rarefied OR not_rarefied information to the sample data
  if (remove_not_rarefied) {
    # Remove not-rarefied samples
    ps_rare <- phyloseq::prune_samples(all_names[rarefy_id], ps_rare)
    # Return rarefied phyloseq object
    message("Not-rarefied samples were removed from output as you specified.")
  } else {
    # Return rarefied phyloseq object
    message2 <- "Rarefied/not-rarefied samples were kept in the phyloseq object."
    message3 <- "Sequence reads of the not-rarefied samples were not changed."
    message4 <- "Please check \'rarefied\' column of sample_data()."
    message(message2); message(message3); message(message4)
  }

  # Return results
  if (include_iNEXT_results) {
    return(list(ps_rare, inext_out))
  } else {
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
plot_rarefy <- function (ps_obj,
                         plot_rarefied_point = TRUE,
                         #se = FALSE,
                         plot_slope = FALSE,
                         linetype = 1) {
  # Extract iNEXT result
  inext_res <- sapply(ps_obj[[2]], `[`, 2)
  names(inext_res) <- phyloseq::sample_names(ps_obj[[1]])

  # Add sample sames to each data
  inext_df <- purrr::map_dfr(inext_res, function(x) return(data.frame(x$size_based)), .id = "sample")
  colnames(inext_df)[which(colnames(inext_df) == "m")] <- "x"
  colnames(inext_df)[which(colnames(inext_df) == "qD")] <- "y"
  #dplyr::rename(x = m, y = qD)

  # Compile rarefied data
  rare_df2 <- data.frame(phyloseq::sample_data(ps_obj[[1]]))
  colnames(rare_df2)[which(colnames(rare_df2) == "rarefied_reads")] <- "x"
  colnames(rare_df2)[which(colnames(rare_df2) == "rarefied_pred_n_taxa")] <- "y"
  #dplyr::rename(x = rarefied_reads, y = rarefied_n_taxa)
  rare_df2$rarefied_slope <- 1 - rare_df2$rarefied_coverage
  rare_df2$sample <- rownames(rare_df2)
  rare_df2_slope <- rare_df2$rarefied_slope
  rare_df2_intercept <- rare_df2$y - rare_df2$rarefied_slope * rare_df2$x

  # Visualize with ggplot
  g_rare <- ggplot2::ggplot(inext_df, ggplot2::aes_(x = as.name("x"), y = as.name("y"), color = as.name("sample"))) +
    ggplot2::geom_line(size = 0.5) +
    ggplot2::xlab("Sequence reads") +
    ggplot2::ylab("The number of species") +
    NULL

  if (plot_rarefied_point) {
    # Add rarefied information
    g_rare <- g_rare +
      ggplot2::geom_point(data = rare_df2,
                          ggplot2::aes_string(x = "x", y = "y"),
                          color = "gray30", size = 3, shape = 18, alpha = 0.6) +
      NULL
  }

  # Add slopes
  if (plot_slope) {
    g_rare <- g_rare +
      ggplot2::geom_abline(slope = rare_df2_slope,
                           intercept = rare_df2_intercept, color = "gray60",
                           linetype = linetype, linewidth = 0.3, alpha = 0.5)
  }

  # Return ggplot object
  return(g_rare)
}


#' @title Estimate the required sample size for a particular coverage
#' @description \code{coverage_to_samplesize} estimates the required sample size for a particular coverage. using iNEXT package
#' @param x Vector of species abundances.
#' @param coverage The desired sample completeness that we want to achieve.
#' @param add_attr Logical. Specify whether attributes are added.
#' @return The required sample size (`mm`).
#' @noRd
#' @keywords internal

# ---------------------------------------------- #
# coverage_to_samplesize
# ---------------------------------------------- #
# From https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_coverage.R#sym-coverage_to_samplesize
#
## Estimate the required sample size for a particular coverage
## the code is based on iNEXT:::invChat.Ind by Johnson Hsieh (d76e3b8, Nov 12, 2016)
# https://github.com/JohnsonHsieh/iNEXT/
coverage_to_samplesize <- function(x, coverage = 0.95, add_attr = F){
  # x = vector of species abundances
  # coverage = the desired sample completeness that we want to achieve
  # iNEXT:::invChat.Ind(x, C = coverage)$m
  ## Total number of reads and the observed sample coverage
  n <- sum(x)
  refC <- Chat.Ind(x, n)

  ## Interpolation
  f <- function(m, C) abs(Chat.Ind(x, m) - C)
  if(refC >= coverage){
    opt <- stats::optimize(f, C = coverage, lower = 0, upper = sum(x))
    mm <- opt$minimum
    mm <- round(mm)
  }

  ## Extrapolation
  if(refC < coverage){
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    if(f1>0 & f2>0)  {A <- (n-1)*f1/((n-1)*f1+2*f2)}
    if(f1>1 & f2==0) {A <- (n-1)*(f1-1)/((n-1)*(f1-1)+2)}
    if(f1==1 & f2==0){A <- 1}
    if(f1==0 & f2==0){A <- 1}
    mm <- (log(n/f1)+log(1-coverage))/log(A)-1
    mm <- n + mm
    mm <- round(mm)
  }

  ## Add attributes
  if(add_attr){
    if(refC > coverage) { attr(mm, "method") <- "interpolated" }
    if(refC <= coverage){ attr(mm, "method") <- "extrapolated" }
    attr(mm, "ObservedCoverage") <- refC
    attr(mm, "RequestedCoverage") <- coverage
  }

  return(mm)
}
