
# macam: A collection of miscellaneous R functions for ecology

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/492736367.svg)](https://zenodo.org/badge/latestdoi/492736367)
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
- `cov_info()`: Get coverage and diversity information of a `phyloseq` object using `iNEXT` package.
- `rarefy_even_coverage()`: Perform coverage-based rarefaction.
- `plot_rarefy()`: Visualize results of coverage-based rarefaction.
- `taxa_name_summarize()`: Summarize taxa names in phyloseq object for visualization purpose.

### Time series analysis
- `s_map_rgl()`: Perform the regularized S-map in [Cenci et al. (2019)](https://doi.org/10.1111/2041-210X.13150). A wrapper function of `extended_lnlp()`.
- `extended_lnlp()`: A generalized function for the regularized S-map. For the multivariate S-map, please use this function.
- `twin_surrogate_cpp()`: Generate twin surrogate time series from an original time series (see [Thiel et al. 2006](https://doi.org/10.1209/epl/i2006-10147-0)).
- `uic_across()`: Perform UIC across columns (for the MDR S-map) (beta version).
- `make_block_mvd()`: Generate data.frame for the calculation of multiview distance (for the MDR S-map) (beta version).
- `compute_mvd()`: Compute multiview distance (for the MDR S-map) (beta version).
- `s_map_mdr()`: Perform the MDR S-map in [Chang et al. (2021)](https://doi.org/10.1111/ele.13897)  (beta version).
- `s_map_mdr_all()`: A convenient all-in-one wrapper function for the MDR S-map. To fine-tune parameters, use step-by-step functions such as `compute_mvd()` etc (beta version).

### ggplot2 functions
- `label_10_to_power()`: Convert a numeric variable to a scientific notation (e.g., 1500 will be $1.5 \times 10^3$)

### data
- `data_4sp`: Demo time series data for time series analysis. Data is from [Ushio et al. (2018)](https://doi.org/10.1038/nature25504).
- `asv_sheet`: A demo ASV sheet for sequence analysis. Data is from [Ushio (2019)](https://doi.org/10.1111/2041-210X.13204).
- `sample_sheet`: A demo sample sheet for sequence analysis. Data is from [Ushio (2019)](https://doi.org/10.1111/2041-210X.13204).
- `tax_sheet`: A demo taxa sheet for sequence analysis. Data is from [Ushio (2019)](https://doi.org/10.1111/2041-210X.13204).


## References
- Cenci et al. (2019) Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.13150
- Chang et al. (2021) Ecology Letters. https://doi.org/10.1111/ele.13897
- Chiu & Chao (2016) PeerJ. https://doi.org/10.7717/peerj.1634
- Hsieh et al. (2016) Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.12613
- Osada & Ushio (2020) rUIC:Unified Information-theoretic Causality for R. https://doi.org/10.5281/zenodo.5163234
- Osada et al. (2023) bioRxiv. https://doi.org/10.1101/2023.04.20.537743
- Mikryukov (2018) vmikk/metagMisc: v.0.0.4. https://doi.org/10.5281/zenodo.571403
