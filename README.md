
# macam: A collection of miscellaneous functions for ecology in R

<!-- badges: start -->
<!-- badges: end -->

This package is a collection of miscellaneous functions for ecological studies.

## License
See LICENSE.


## Installation

You can install the development version of macam from [GitHub](https://github.com/) with:

``` r
#install.packages("remotes")
remotes::install_github("ong8181/macam")

# Or build_vignettes = TRUE to install with vignettes
remotes::install_github("ong8181/macam", build_vignettes = TRUE, force = TRUE)
```

## Functions
Type `?macam::xxxxx()` for detail. Also please check `browseVignette("macam")`.

### File/directory handling
- `outdir_create()`: Create output directory.
- `save_workspace()`: Save workspace.
- `save_session_info()`: Save session information.

### Sequence analysis
- `AllOrients()`: List up all orientations of primer sequence (from [DADA2 tutrial](https://benjjneb.github.io/dada2/ITS_workflow.html)).
- `PrimerHits()`: Count the number of primer hits (from [DADA2 tutrial](https://benjjneb.github.io/dada2/ITS_workflow.html))
- `rarefy_even_coverage()`: Perform coverage-based rarefaction.
- `plot_rarefy()`: Visualize results of coverage-based rarefaction.
- `taxa_name_summarize()`: Summarize taxa names in phyloseq object for visualization purpose.

### Time series analysis
- `s_map_rgl()`: Perform the regularized S-map in [Cenci et al. (2019)](https://doi.org/10.1111/2041-210X.13150). A wrapper function of `extended_lnlp()`.
- `extended_lnlp()`: A generalized function for the regularized S-map. For the multivariate S-map, please use this function.
- `twin_surrogate_cpp()`: Generate twin surrogate time series from an original time series (see [Thiel et al. 2006](https://doi.org/10.1209/epl/i2006-10147-0)).
- `uic_across()`: Perform UIC across columns.
- `make_block_mvd()`: Generate data.frame for the calculation of multiview distance
- `compute_mvd()`: Compute multiview distance
- `s_map_mdr()`: Perform the MDR S-map in [Chang et al. (2021)](https://doi.org/10.1111/ele.13897).
- `s_map_mdr_all()`: A convenient all-in-one wrapper function for the MDR S-map. To fine-tune parameters, use step-by-step functions such as `compute_mvd()` etc.

### ggplot2 functions
- `label_10_to_power()`: Convert a numeric variable to a scientific notation (e.g., 1500 will be $1.5 \times 10^3$)
