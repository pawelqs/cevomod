
#' Get sample mutation burden
#' @export
get_sample_mutation_burden <- function(cd) {
  SNVs(cd) |>
    filter(.data$VAF > 0, !is.na(.data$VAF)) |>
    group_by(.data$sample_id) |>
    summarise(mutation_burden = n(), .groups = "drop")
}


#' Get tumor mutation burden
#' @export
get_patient_mutation_burden <- function(cd) {
  SNVs(cd) |>
    filter(.data$VAF > 0, !is.na(.data$VAF)) |>
    group_by(.data$patient_id) |>
    summarise(mutation_burden = n(), .groups = "drop")
}

