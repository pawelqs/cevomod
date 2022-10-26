
get_sequencing_depths <- function(object, ...) {
  snvs <- SNVs(object)
  sequencing_depths <- snvs |>
    mutate(VAF = round(.data$VAF, digits = 2)) |>
    group_by(.data$sample_id, .data$VAF) |>
    summarise(
      mean_DP = mean(.data$DP, na.rm = TRUE),
      median_DP = stats::median(.data$DP, na.rm = TRUE),
      sd_DP = sd(.data$DP, na.rm = TRUE),
      .groups = "drop"
    )
  sequencing_depths
}
