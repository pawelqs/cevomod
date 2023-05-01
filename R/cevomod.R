
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

#' TCGA BRCA cevodata dataset
#' @name tcga_brca
#' @docType data
NULL

#' Small TCGA BRCA cevodata dataset
#' @name tcga_brca_test
#' @docType data
NULL

#' 4 test samples
#' @name test_data
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
  object <- object |>
    calc_mutation_frequencies() |>
    calc_SFS() |>
    calc_cumulative_tails() |>
    calc_Mf_1f() |>
    fit_williams_neutral_models() |>
    fit_subclones() |>
    fit_powerlaw_tail_optim() |>
    fit_subclones()
  object$active_models <- "williams_neutral_subclones"
  object
}


