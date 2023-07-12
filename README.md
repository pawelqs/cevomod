
# cevomod <img src="man/figures/logo.png" align="right" height="120" />

<!-- badges: start -->
<!--[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) -->
[![R-CMD-check](https://github.com/pawelqs/cevomod/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pawelqs/cevomod/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/pawelqs/cevomod/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/pawelqs/cevomod/actions/workflows/test-coverage.yaml)
<!-- badges: end -->


*cevomod* is a package that implements methods for analyzing cancer evolution from the Next Generation Sequencing data.

The modeling approach implemented in *cevomod* was inspired by the [MOBSTER](https://caravagnalab.github.io/mobster/index.html) package, which models the distribution of Variant Allele Frequencies in the sample with a mixture of power-law-shaped and binomial distributions. However, MOBSTER fails to recognize the power-law component (the so-called neutral tail) in Whole Exome Sequencing data or data with insufficient sequencing coverage. cevomod implements methods that can fit the model to the data with significant loss of neutral tail variants.


![](man/figures/introduction_figure.svg)


## Installation

You can install the development version of cevomod from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pawelqs/cevomod")
```

## Help and support

[GitHub Issues](https://github.com/pawelqs/cevomod/issues)

