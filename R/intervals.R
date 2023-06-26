
#' Get mutation frequency intervals
#'
#' Cuts requested column e.g. CCF or VAF into intervals. Adds f_interval and f
#' columns to the input SNVs tibble
#'
#' @param snvs SNVs tibble
#' @param bins Number of bins
#' @param column Name of column with frequencies, eg. VAF or CCF
#'
#' @return Tibble with f_interval and f columns
#' @export
cut_f_intervals <- function(snvs, column, bins = NULL) {
  breaks <- get_interval_breaks(snvs, bins = bins) |>
    enframe(name = "sample_id", value = "breaks")
  res <- snvs |>
    nest_by(.data$sample_id) |>
    left_join(breaks, by = "sample_id") |>
    mutate(data = list(cut_f(.data$data, .data$breaks, column = column))) |>
    ungroup()
  interval_levels <- res |>
    select("sample_id", "data") |>
    deframe() |>
    map("f_interval") |>
    map(levels)
  res$data <- res$data |>
    map(~mutate(.x, across(f_interval, as.character)))
  res <- res |>
    select("sample_id", "data") |>
    unnest("data") |>
    mutate(f = get_interval_centers(.data$f_interval), .after = "f_interval")
  attr(res, "intervals") <- interval_levels
  res
}


cut_f <- function(tbl, breaks, column) {
  tbl |>
    mutate(
      f_interval = cut(.data[[column]], breaks = breaks),
      .after = all_of(column)
    )
}


get_interval_breaks <- function(object, bins = NULL, sample_id = NULL) {
  if (is.null(bins)) {
    bins_by_sample <- get_sample_sequencing_depths(object) |>
      transmute(
        .data$sample_id,
        bins = round(.data$median_DP),
        bins = if_else(.data$bins > 100, 100, bins)
      )
  } else {
    bins_by_sample <- tibble(
      sample_id = unique(object$sample_id),
      bins = bins
    )
  }
  breaks <- bins_by_sample |>
    deframe() |>
    map(~c(-1/.x, seq(0, 1, length.out = .x + 1)))

  if (is.null(sample_id)) {
    breaks
  } else {
    breaks[[sample_id]]
  }
}


complete_missing_f_intervals <- function(tbl, intervals) {
  intervals <- intervals |>
    enframe(name = "sample_id", value = "f_interval") |>
    unnest("f_interval")
  missing_intervals <- intervals |>
    anti_join(tbl, by = c("sample_id", "f_interval"))
  tbl |>
    bind_rows(missing_intervals) |>
    arrange(.data$sample_id, .data$f_interval)
}


get_interval_centers <- function(intervals) {
  res <- intervals |>
    str_replace_all("[\\(\\)\\[\\]]", "") |>
    as_tibble_col("from_and_to") |>
    separate_wider_delim("from_and_to", names = c("from", "to"), delim = ",") |>
    map_df(parse_double) |>
    mutate(centers = .data$from + (.data$to - .data$from) / 2)
  res$centers
}


get_interval_width <- function(intervals) {
  res <- intervals |>
    str_replace_all("[\\(\\)\\[\\]]", "") |>
    as_tibble_col("from_and_to") |>
    separate_wider_delim("from_and_to", names = c("from", "to"), delim = ",") |>
    map_df(parse_double) |>
    mutate(width = .data$to - .data$from)
  stats::median(res$width)
}

