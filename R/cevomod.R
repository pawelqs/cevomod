
#' Tbl of driver genes from Bailey, Ding et al. 'Comprehensive Characterization
#' of Cancer Driver Genes and Mutations', Cell, 2-18
#' https://doi.org/10.1016/j.cell.2018.02.060
#' @name driver_genes
#' @docType data
NULL

#' Custom Variant Classification
#' @name variant_classification
#' @docType data
NULL

#' TCGA-BRCA SNVs dataset
#' @name snvs_tcga_brca
#' @docType data
NULL

#' Small TCGA BRCA cevodata dataset
#' @name tcga_brca_test
#' @docType data
NULL


#' Run cevodata pipeline
#' @param object cevodata object
#' @param ... other args
#' @export
run_cevomod <- function(object, ...) {
  UseMethod("run_cevomod")
}


#' @export
run_cevomod.cevodata <- function(object, ...) {
  object |>
    calc_mutation_frequencies() |>
    calc_SFS() |>
    calc_cumulative_tails() |>
    calc_Mf_1f() |>
    fit_williams_neutral_models() |>
    fit_subclones() |>
    fit_tung_durrett_models() |>
    fit_subclones()
  object$active_models <- "williams_neutral_subclones"
  object
}


