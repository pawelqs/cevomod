
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
#' tcga_brca_test
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
  snvs <- cut_VAF_intervals(object, bins = bins)
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


#' Plot SFS
#'
#' @param x tibble with calc_SFS() results
#' @param ... futher passed to geom_()
#' @param geom geom
#' @return ggplot obj
#' @export
plot.cevo_SFS_tbl <- function(x, ..., geom = "bar") {
  if (geom == "bar") {
    x <- x |>
      group_by(.data$sample_id) |>
      mutate(width = 0.9 / n())
    geom <- list(
      geom_bar(
        aes(fill = .data$sample_id, width = .data$width),
        stat = "identity", alpha = 0.8, ...
      ),
      facet_wrap(~.data$sample_id, scales = "free")
    )
  } else if (geom == "line") {
    geom <- list(
      geom_line(...)
    )
  }

  x |>
    filter(.data$VAF >= 0) |>
    ggplot(aes(.data$VAF, .data$y, color = .data$sample_id)) +
    geom +
    theme_ellie(n = n_distinct(x$sample_id)) +
    hide_legend() +
    labs(
      title = "SFS",
      y = "count"
    )
}


#' @describeIn sfs Plot SFS
#' @export
plot_SFS.cevodata <- function(object, ..., geom = "bar") {
  plot(object$models$SFS, geom = geom)
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
