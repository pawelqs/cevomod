
## cevomod 2.3.0
* subclones can be fitted using [CliP](https://github.com/wwylab/CliP) [(Jiang et al., CliP: subclonal architecture reconstruction of cancer cells in DNA sequencing data using a penalized likelihood model, bioRxiv 2021)](https://www.biorxiv.org/content/10.1101/2021.03.31.437383v1) (`fit_subclones(method = "CliP")`). This option requires that the Apptainer is installed. CliP container image can be build with `build_clip_container()`

## cevomod 2.2.0
* fit_powerlaw_tail_fixed() has a bootstrap option

## cevomod 2.1.0
* cevomod is integrated with a helper [readthis](https://pawelqs.github.io/readthis/index.html) package, designed for bulk reading of variant files from algorithms such as Mutect2, Strelka, ASCAT, or FACETS, in the cevomod-friendly data format. Objects returned by `readthis::read_*()` functions can be added to the cevodata object using a general `add_data()` function.


## cevomod 2.0.0
* cevomod functions can no utilize VAF or CCF (Cancer Cell Fraction) as a measure
  of mutation frequency. CCF is calculated using the formula introduced in [Dentro et al. *Principles of Reconstructing the Subclonal Architecture of Cancers* (2015)](https://doi.org/10.1101/cshperspect.a026625)


## cevomod 1.1.0
* cevodata export to [CliP](https://github.com/wwylab/CliP) implemented


## cevomod 1.0.0
* cevodata class implementation
* fitting the power-law tails with exponent equal to 2 using $M(f) \sim 1/f$ statistic
* fitting the power-law tails with optimized exponent
* fitting subclones using mclust
* fitting subclones using BMix
* calculation of the evolutionary parameters using the [Williams et al. (2018)](https://doi.org/10.1038/s41588-018-0128-6) equations and the [MOBSTER code](https://github.com/caravagnalab/mobster/blob/master/R/evodynamics.R) [(Caravagna et al. (2020))](https://doi.org/10.1038/s41588-020-0675-5)
