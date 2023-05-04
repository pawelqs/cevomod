


#' Prepare SNVs for analyses
#' @param object cevodata obj
#' @param which which SNVs to use
#' @param bins number of VAF interval bins
#' @param verbose verbose?
#' @export
prepare_SNVs <- function(object, which = default_SNVs(object), bins = NULL, verbose = TRUE) {
  snvs <- SNVs(object, which) |>
    cut_f_intervals(bins = bins)
  object |>
    add_SNV_data(snvs, name = which) |>
    calc_SFS()
}



#' Calc mutation frequencies
#'
#' Currently no method is implemented and this function only initializes
#' f column in active SNVs slot with values from VAF column
#'
#' @name mutation_frequencies

#' @rdname mutation_frequencies
#' @export
calc_mutation_frequencies <- function(object, ...) {
  UseMethod("calc_mutation_frequencies")
}


#' @rdname mutation_frequencies
#' @param object data object
#' @param method method: "use_VAF". Other methods will be implemented later
#' @param ... other arguments
#' @export
calc_mutation_frequencies.cevodata <- function(object, method = "use_VAF", ...) {
  if (method == "use_VAF") {
    snvs <- SNVs(object) |>
      mutate(f = .data$VAF, .after = "VAF")
  }
  object |>
    add_SNV_data(snvs, name = object$active_SNVs)
}
