
#' Get sample mutation burden
#' @param object object
#' @param snvs which SNVs to use
#' @param ... other arguments (currently nut used)
#' @export
get_sample_mutation_burden <- function(object, ...) {
  UseMethod("get_sample_mutation_burden")
}


#' @rdname get_sample_mutation_burden
#' @export
get_sample_mutation_burden.cevodata <- function(object, snvs = default_SNVs(object), ...) {
  SNVs(object, which = snvs) |>
    get_sample_mutation_burden()
}


#' @rdname get_sample_mutation_burden
#' @export
get_sample_mutation_burden.cevo_snvs <- function(object, ...) {
  object |>
    filter(.data$VAF > 0, !is.na(.data$VAF)) |>
    group_by(.data$sample_id) |>
    summarise(mutation_burden = n(), .groups = "drop")
}


#' Get tumor mutation burden
#' @param object cevodata object
#' @param snvs snvs slot
#' @export
get_patient_mutation_burden <- function(object, snvs = default_SNVs(object)) {
  SNVs(object, which = which) |>
    filter(.data$VAF > 0, !is.na(.data$VAF)) |>
    left_join(object$metadata, by = "sample_id") |>
    group_by(.data$patient_id) |>
    summarise(mutation_burden = n(), .groups = "drop")
}

