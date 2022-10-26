
get_sequencing_depths <- function(object, ...) {
  snvs <- SNVs(object)
  sequencing_depths <- snvs |>
    mutate(VAF = round(VAF, digits = 2)) |>
    group_by(sample_id, VAF) |>
    summarise(
      mean_DP = mean(DP, na.rm = TRUE),
      median_DP = median(DP, na.rm = TRUE),
      sd_DP = sd(DP, na.rm = TRUE),
      .groups = "drop"
    )
  sequencing_depths
}
