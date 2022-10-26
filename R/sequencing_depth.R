
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

  completed_depths <- sequencing_depths |>
    nest_by(.data$sample_id) |>
    summarise(
      VAF = 1:100/100,
      mean_DP = stats::approx(.data$data$VAF, .data$data$mean_DP, xout = .data$VAF, rule = 2)$y,
      median_DP = stats::approx(.data$data$VAF, .data$data$median_DP, xout = .data$VAF, rule = 2)$y,
      sd_DP = stats::approx(.data$data$VAF, .data$data$sd_DP, xout = .data$VAF, rule = 2)$y,
      .groups = "drop"
    )

  completed_depths
}


