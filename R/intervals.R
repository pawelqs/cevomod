
#' Cut VAF into intervals
#' @param object object
#' @param ... other params
#' @export
cut_VAF_intervals <- function(object, ...) {
  UseMethod("cut_VAF_intervals")
}


#' @rdname cut_VAF_intervals
#' @param bins number of bins
#' @export
cut_VAF_intervals.cevo_snvs <- function(object, bins = NULL, ...) {
  breaks <- get_interval_breaks(object, bins = bins) |>
    enframe(name = "sample_id", value = "breaks")
  res <- object |>
    nest_by(.data$sample_id) |>
    left_join(breaks, by = "sample_id") |>
    mutate(data = list(cut_VAF(.data$data, .data$breaks))) |>
    ungroup()
  interval_levels <- res |>
    select("sample_id", "data") |>
    deframe() |>
    map("VAF_interval") |>
    map(levels)
  res <- res |>
    select("sample_id", "data") |>
    unnest("data")
  attr(res, "intervals") <- interval_levels
  res
}


cut_VAF <- function(tbl, breaks) {
  tbl |>
    mutate(
      VAF_interval = cut(.data$VAF, breaks = breaks),
      .before = "VAF"
    )
}


get_interval_breaks <- function(object, bins = NULL) {
  if (is.null(bins)) {
    bins_by_sample <- get_sample_sequencing_depths(object) |>
      transmute(.data$sample_id, bins = round(.data$median_DP))
  } else {
    bins_by_sample <- tibble(
      sample_id = unique(object$sample_id),
      bins = bins
    )
  }
  breaks <- bins_by_sample |>
    deframe() |>
    map(~c(-1/.x, seq(0, 1, length.out = .x + 1)))
  breaks
}


complete_missing_VAF_intervals <- function(tbl, intervals) {
  intervals <- intervals |>
    enframe(name = "sample_id", value = "VAF_interval") |>
    unnest("VAF_interval")
  missing_intervals <- intervals |>
    anti_join(tbl, by = c("sample_id", "VAF_interval"))
  tbl |>
    bind_rows(missing_intervals) |>
    arrange(.data$sample_id, .data$VAF_interval)
}


get_interval_centers <- function(intervals) {
  intervals |>
    str_replace_all("[\\(\\)\\[\\]]", "") |>
    str_split(pattern = ",") |>
    map(parse_double) |>
    map(mean) |>
    unlist()
}


get_interval_width <- function(intervals) {
  intervals |>
    str_replace_all("[\\(\\)\\[\\]]", "") |>
    str_split(pattern = ",") |>
    map(parse_double) |>
    map(~.x[[2]] - .x[[1]]) |>
    unlist() |>
    stats::median()
}

