

#' Small TCGA BRCA cevodata dataset
#' @name tcga_brca_fitted
#' @docType data
NULL

#' 4 test samples
#' @name test_data
#' @docType data
NULL


#' cevomod results for 4 test samples
#' @name test_data_fitted
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
    fit_powerlaw_tail_fixed() |>
    fit_subclones() |>
    fit_powerlaw_tail_optim() |>
    fit_subclones()
  active_models(object) <- "williams_neutral_subclones"
  object
}


