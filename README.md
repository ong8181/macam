
# macam: A collection of miscellaneous functions for ecology in R

<!-- badges: start -->
<!-- badges: end -->

This pacakge is a collection of miscellaneous functions for ecological studies.

## Installation

You can install the development version of macam from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ong8181/macam")
```

## Functions
### File/directory handling
- `outdir_create()`: Create output directory.
- `save_workspace()`: Save workspace.
- `save_session_info()`: Save session information.

### Sequence analysis
- `AllOrients()`: List up all orientations of primer sequence (from [DADA2 tutrial](https://benjjneb.github.io/dada2/ITS_workflow.html)).
- `PrimerHits()`: Count the number of primer hits (from [DADA2 tutrial](https://benjjneb.github.io/dada2/ITS_workflow.html))
- `rarefy_even_coverage()`: Perform coverage-based rarefaction).
- `plot_rarefy()`: Visualize results of coverage-based rarefaction.
