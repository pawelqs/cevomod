
#' Site Frequency Spectra
#'
#' Creates  cevodata$models$SFS with the groupping variables and:
#'   - n columnt with the number of mutations in the VAF interval
#'   - x and y columns describing SFS
#'   - y_scaled with y values scaled to the range 0-1
#'
#' @param object SNVs tibble object
#' @param bins resolution of the cumulative tails calculation
#' @param geom geom
#' @param ... other arguments
#' @examples
#' data("tcga_brca_test")
#'
#' tcga_brca_test |>
#'   calc_SFS() |>
#'   plot_SFS() +
#'   layer_mutations(tcga_brca_test, drivers = "BRCA")
#' @name sfs
NULL


#' @describeIn sfs Calculate SFS
#' @export
calc_SFS <- function(object, ...) {
  UseMethod("calc_SFS")
}


#' @describeIn sfs Calculate SFS
#' @export
calc_SFS.cevodata <- function(object, bins = NULL, ...) {
  object$models[["SFS"]] <- SNVs(object) |>
    calc_SFS(bins = bins)
  object
}


#' @describeIn sfs Calculate SFS
#' @export
calc_SFS.cevo_snvs <- function(object, bins = NULL, ...) {
  snvs <- cut_f_intervals(object, bins = bins) |>
    filter(.data$VAF > 0)
  intervals <- attributes(snvs)$intervals
  res <- snvs |>
    group_by(.data$sample_id, .data$VAF_interval) |>
    summarise(y = n(), .groups = "drop_last") |>
    complete_missing_VAF_intervals(intervals) |>
    replace_na(list(y = 0)) |>
    mutate(VAF = get_interval_centers(.data$VAF_interval), .after = "VAF_interval") |>
    mutate(y_scaled = round(.data$y / sum(.data$y), digits = 4)) |>
    ungroup()
  class(res) <- c("cevo_SFS_tbl", class(res))
  res
}


#' @describeIn sfs Plot SFS
#' @export
plot_SFS <- function(object, ...) {
  UseMethod("plot_SFS")
}


#' @describeIn sfs Plot SFS
#' @param mapping aes()
#' @export
plot_SFS.cevodata <- function(object, mapping = NULL, ..., geom = "bar") {
  if (is.null(object$models$SFS)) {
    object <- calc_SFS(object)
  }
  # TODO: Fix 'width' warning
  object$models$SFS |>
    left_join(object$metadata, by = "sample_id") |>
    plot(mapping = mapping, ..., geom = geom)
}


#' Plot SFS
#'
#' @param x tibble with calc_SFS() results
#' @param mapping aes()
#' @param ... futher passed to geom_()
#' @param geom geom
#' @return ggplot obj
#' @export
plot.cevo_SFS_tbl <- function(x, mapping = NULL, ..., geom = "bar") {
  default_mapping <- aes(.data$VAF, .data$y, group = .data$sample_id)

  if (geom == "bar") {
    x <- x |>
      group_by(.data$sample_id) |>
      mutate(width = 0.9 / n())
    bar_mapping <- aes(width = .data$width)
    p <- ggplot(x) +
      join_aes(default_mapping, mapping) +
      geom_bar(
        join_aes(bar_mapping, mapping),
        stat = "identity", alpha = 0.8, ...
      ) +
      facet_wrap(~.data$sample_id, scales = "free")
  } else if (geom == "line") {
    p <- ggplot(x, join_aes(default_mapping, mapping)) +
      geom_line(...)
  }

  p + labs(title = "SFS", y = "count")
}


#' SFS stat
#'
#' Modification of stat_bin with custom parameters
#' @param mapping aes()
#' @param data data
#' @param binwidth binwidth
#' @param geom geom
#' @param position position
#' @param ... args passed to stat_bin
#'
#' @export
stat_SFS <- function(mapping = NULL, data = NULL,
                     binwidth = 0.01, geom = "line", position = "identity", ...) {
  stat_bin(mapping = mapping, data = data, binwidth = binwidth, geom = geom, position = position, ...)
}


#' @describeIn sfs Get SFS
#' @export
get_SFS <- function(object, ...) {
  UseMethod("get_SFS")
}


#' @describeIn sfs Get SFS
#' @param model_name name of slot with SFS statistics
#' @param verbose verbose?
#' @export
get_SFS.cevodata <- function(object, model_name = "SFS", verbose = TRUE, ...) {
  sfs <- object$models[[model_name]]
  if (is.null(sfs)) {
    msg(
      "SFS's not calculated yet. Calculating with sample DP as number of bins",
      verbose = verbose
    )
    object <- calc_SFS(object)
    sfs <- object$models[[model_name]]
  }
  sfs
}


get_non_zero_SFS_range <- function(sfs,
                                   allowed_zero_bins = 1,
                                   y_treshold = 1,
                                   y_threshold_pct = 0.01) {
  sfs |>
    group_by(.data$sample_id) |>
    mutate(
      empty_bin = (.data$y < y_treshold) | (.data$y < max(.data$y) * y_threshold_pct),
      segment_number = segment(.data$empty_bin),
      empty_low_VAF_range = .data$empty_bin & .data$segment_number == 0
    ) |>
    filter(!.data$empty_low_VAF_range) |>
    group_by(.data$sample_id, .data$segment_number) |>
    mutate(
      segment_length = n(),
      keep = !.data$empty_bin | (.data$empty_bin & (.data$segment_length <= allowed_zero_bins))
    ) |>
    group_by(.data$sample_id) |>
    mutate(new_segments = segment(.data$keep)) |>
    filter(.data$new_segments == 0) |>
    summarise(
      from = min(.data$VAF),
      to = max(.data$VAF)
    )
}
