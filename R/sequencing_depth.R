
#' Get sample sequencing depths
#' @param object object
#' @param ... other args
#' @export
get_sample_sequencing_depths <- function(object, ...) {
  UseMethod("get_sample_sequencing_depths")
}


#' @export
get_sample_sequencing_depths.cevodata <- function(object, ...) {
  SNVs(object) |>
    get_sample_sequencing_depths()
}


#' @export
get_sample_sequencing_depths.cevo_snvs <- function(object, ...) {
  sequencing_depths <- object |>
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
    mutate(f = round(.data$f, digits = 2)) |>
    group_by(.data$sample_id, .data$f) |>
    summarise(
      mean_DP = mean(.data$DP, na.rm = TRUE),
      median_DP = stats::median(.data$DP, na.rm = TRUE),
      sd_DP = sd(.data$DP, na.rm = TRUE),
      .groups = "drop"
    )

  completed_depths <- sequencing_depths |>
    nest_by(.data$sample_id) |>
    reframe(
      f = 1:100/100,
      mean_DP = stats::approx(.data$data$f, .data$data$mean_DP, xout = .data$f, rule = 2)$y,
      median_DP = stats::approx(.data$data$f, .data$data$median_DP, xout = .data$f, rule = 2)$y,
      sd_DP = stats::approx(.data$data$f, .data$data$sd_DP, xout = .data$f, rule = 2)$y
    )

  completed_depths
}


