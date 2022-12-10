
get_SNVs_wider <- function(object, fill_na = NULL) {
  patients_to_samples <- object$metadata |>
    select(patient_id:sample)

  snvs <- SNVs(object) |>
    select(sample_id, chrom:alt, VAF) |>
    left_join(patients_to_samples, by = "sample_id") |>
    unite(mutation_id, chrom:alt, sep = "-") |>
    select(-sample_id, -patient_id) |>
    pivot_wider(names_from = "sample", values_from = "VAF")

  if (!is.null(fill_na)) {
    snvs[is.na(snvs)] <- 0
  }
  snvs
}
