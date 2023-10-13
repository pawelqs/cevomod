
#' Site Frequency Spectra
#'
#' Site Frequency Spectra (or Variant Allele Frequency Spectra) are the main
#' statistic used by cevomod.
#'
#' @param object SNVs tibble object
#' @param which_snvs Which SNVs to use?
#' @param column VAF or CCF/2
#' @param bins Resolution of the cumulative tails calculation
#' @param verbose Verbose?
#' @param geom Geom
#' @param ... other arguments
#' @examples
#' data("test_data")
#'
#' test_data |>
#'   calc_SFS() |>
#'   plot_SFS() +
#'   layer_mutations(test_data, mapping = ggplot2::aes(x = VAF), drivers = "BRCA")
#' @name sfs
NULL


#' @describeIn sfs Calculates spectra for all samples and saves and saves them
#' in cevodata$models$SFS tibble.
#'
#' SFS columns description:
#' - y number of mutations in the frequency interval
#' - y_scaled with y values scaled to the range 0-1
#'
#' @export
calc_SFS <- function(object, ...) {
  UseMethod("calc_SFS")
}


#' @describeIn sfs method for cevodata object
#' @export
calc_SFS.cevodata <- function(object,
                              which_snvs = default_SNVs(object),
                              column = get_frequency_measure_name(object, which_snvs),
                              bins = NULL,
                              verbose = get_cevomod_verbosity(),
                              ...) {
  object$models[["SFS"]] <- SNVs(object, which_snvs) |>
    calc_SFS(column = column, bins = bins, verbose = verbose)
  object
}


#' @describeIn sfs method for cevo_snvs object
#' @export
calc_SFS.cevo_snvs <- function(object,
                               column = get_frequency_measure_name(object),
                               bins = NULL,
                               verbose = get_cevomod_verbosity(),
                               ...) {
  msg("Calculating SFS statistics", verbose = verbose)

  if (is.null(object[["f_interval"]]) | !is.null(bins)) {
    msg("Calculating f intervals, using ", column, " column", verbose = verbose)
    snvs <- cut_f_intervals(object, column = column, bins = bins) |>
      filter(.data$f > 0)
  } else {
    snvs <- object |>
      filter(.data$f > 0)
  }
  intervals <- attributes(snvs)$intervals
  res <- snvs |>
    group_by(.data$sample_id, .data$f_interval) |>
    summarise(y = n(), .groups = "drop_last") |>
    complete_missing_f_intervals(intervals) |>
    replace_na(list(y = 0)) |>
    mutate(f = get_interval_centers(.data$f_interval), .after = "f_interval") |>
    mutate(y_scaled = round(.data$y / sum(.data$y), digits = 4)) |>
    ungroup()
  class(res) <- c("cevo_SFS_tbl", class(res))
  attr(res, "f_column") <- attributes(snvs)[["f_column"]]
  res
}


calc_SFS_resamples <- function(cd, times, verbose = get_cevomod_verbosity()) {
  rlang::check_installed("rsample", reason = "to perform bootstrap sampling of SNVs")
  msg("Splitting SNVs by sample_id", verbose = verbose)

  splitted_snvs <- SNVs(cd) |>
    nest_by(.data$sample_id, .keep = TRUE) |>
    deframe() |>
    map(as_cevo_snvs)

  msg("Resampling SNVs and calculating SFSs", verbose = verbose)
  pb <- if (verbose) progress_bar$new(total = length(splitted_snvs)) else NULL
  pass_verbose <- verbose_down(verbose)

  resamples <- splitted_snvs |>
    map(function(snvs) {
      if (!is.null(pb)) pb$tick()
      resamples <- rsample::bootstraps(snvs, times = times)
      resamples$sfs <- resamples$splits |>
        map(rsample::analysis) |>
        map(intervalize_mutation_frequencies, verbose = pass_verbose) |>
        map(calc_SFS, verbose = pass_verbose)
      resamples$splits <- NULL
      class(resamples) <- c("cevo_SFS_bootstraps", class(resamples))
      resamples
    })

  resamples
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
    mutate(sample_id = parse_factor(.data$sample_id, levels = object$metadata$sample_id)) |>
    plot(mapping = mapping, ..., geom = geom)
}


#' Plot SFS
#'
#' @param x tibble with calc_SFS() results
#' @param mapping aes()
#' @param alpha alpha
#' @param ... futher passed to geom_()
#' @param geom geom
#' @return ggplot obj
#' @export
plot.cevo_SFS_tbl <- function(x, mapping = NULL, alpha = 0.8, ..., geom = "bar") {
  default_mapping <- aes(.data$f, .data$y, group = .data$sample_id)
  x_label <- if (is.null(attr(x, "f_column"))) "f" else attr(x, "f_column")

  if (geom == "bar") {
    x <- x |>
      group_by(.data$sample_id) |>
      mutate(width = 0.9 / n())
    bar_mapping <- aes(width = .data$width)
    p <- ggplot(x) +
      join_aes(default_mapping, mapping) +
      geom_bar(
        join_aes(bar_mapping, mapping),
        stat = "identity", alpha = alpha, ...
      ) +
      facet_wrap(~.data$sample_id, scales = "free")
  } else if (geom == "line") {
    p <- ggplot(x, join_aes(default_mapping, mapping)) +
      geom_line(...)
  }

  p + labs(title = "SFS", y = "count", x = x_label)
}


#' @describeIn sfs Get SFS
#' @param model_name name of slot with SFS statistics
#' @param verbose verbose?
#' @export
get_SFS <- function(object, model_name = "SFS", verbose = TRUE, ...) {
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


#' Get range of non empty SFS bins
#' @param sfs SFS
#' @param allowed_zero_bins number of allowed empty bins in the interval
#' @param y_treshold bins with less mutations will be considered empty
#' @param y_threshold_pct bins that have less mutations than this param times the
#'   height of the higherst peak will be considered empty
#' @keywords internal
get_non_zero_SFS_range <- function(sfs,
                                   allowed_zero_bins = 1,
                                   y_treshold = 1,
                                   y_threshold_pct = 0.01) {
  sfs |>
    group_by(.data$sample_id) |>
    mutate(
      empty_bin = (.data$y < y_treshold) | (.data$y < max(.data$y) * y_threshold_pct),
      segment_number = segment(.data$empty_bin),
      empty_low_f_range = .data$empty_bin & .data$segment_number == 0
    ) |>
    filter(!.data$empty_low_f_range) |>
    group_by(.data$sample_id, .data$segment_number) |>
    mutate(
      segment_length = n(),
      keep = !.data$empty_bin | (.data$empty_bin & (.data$segment_length <= allowed_zero_bins))
    ) |>
    group_by(.data$sample_id) |>
    mutate(new_segments = segment(.data$keep)) |>
    filter(.data$new_segments == 0) |>
    summarise(
      from = min(.data$f),
      to = max(.data$f)
    )
}
