
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


## Last changes
* **v2.1.0** - cevomod is integrated with a helper [readthis](https://pawelqs.github.io/readthis/index.html) package, designed for bulk reading of variant files from algorithms such as Mutect2, Strelka, ASCAT, or FACETS, in the cevomod-friendly data format. Objects returned by `readthis::read_*()` functions can be added to the cevodata object using a general `add_data()` function.
* **v2.0.0** - Starting with version 2.0.0, cevomod can use either VAF or CCF (Cancer Cell Fraction) as a measure of mutation frequency. CCF is a measure of mutation frequency corrected for tumor purity and copy number alterations. CCF can be calculated prior to mutation frequency intervalization using the `calc_mutation_frequencies()` function and requires information on total copy number in tumor and normal tissue and sample purity (tumor cell content). See the Vignettes for more examples.

To see the previous changes in the package see the [Changelog](https://pawelqs.github.io/cevomod/news/index.html)


## Help and support

[GitHub Issues](https://github.com/pawelqs/cevomod/issues)

