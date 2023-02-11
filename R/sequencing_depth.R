
get_sample_sequencing_depths <- function(snvs, ...) {
  sequencing_depths <- snvs |>
    group_by(.data$sample_id) |>
    summarise(
      mean_DP = mean(.data$DP, na.rm = TRUE),
      median_DP = stats::median(.data$DP, na.rm = TRUE),
      sd_DP = sd(.data$DP, na.rm = TRUE),
      .groups = "drop"
    )
  sequencing_depths
}


get_local_sequencing_depths <- function(snvs, ...) {
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
    reframe(
      VAF = 1:100/100,
      mean_DP = stats::approx(.data$data$VAF, .data$data$mean_DP, xout = .data$VAF, rule = 2)$y,
      median_DP = stats::approx(.data$data$VAF, .data$data$median_DP, xout = .data$VAF, rule = 2)$y,
      sd_DP = stats::approx(.data$data$VAF, .data$data$sd_DP, xout = .data$VAF, rule = 2)$y
    )

  completed_depths
}


