
#' Annotate chromosome ploidies in CNV data
#'
#' Adds the normal_cn column to the data. This column is required for e.g.
#' by Dentro CCF calculation method. Requires 'sex' column in the metadata.
#' Males should be encoded by 'M' or "male'.
#'
#' @param object <cevodata> object
#' @param which_cnvs Name of the CNVs slot
#'
#' @return <cevodata> object
#' @export
annotate_normal_cn <- function(object, which_cnvs = default_CNVs(object)) {
  msg("Assuming human genome")
  cnvs <- CNVs(object, which_cnvs) |>
    left_join(get_patient_sex(object), by = "sample_id") |>
    mutate(
      normal_cn = if_else(
        .data$sex %in% c("male", "M") & .data$chrom %in% c("chrX", "chrY"), 1, 2
      )
    )
  object |>
    add_CNV_data(cnvs, name = which_cnvs)
}
