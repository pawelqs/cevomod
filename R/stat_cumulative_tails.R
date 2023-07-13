
#' Cumulative tails
#'
#' cumulative_tails columns:
#'   - f
#'   - n column with the number of mutations in the f interval
#'   - y cumulative tail value
#'   - y_scaled with y values scaled to the range 0-1
#'
#' @param object SNVs tibble object
#' @param which_snvs Which SNVs to use?
#' @param column VAF or CCF/2
#' @param bins Resolution of the cumulative tails calculation
#' @param verbose Verbose?
#' @param ... other arguments
#' @examples
#' data("test_data")
#' test_data |>
#'   calc_cumulative_tails()
#'
#' test_data |>
#'   plot_cumulative_tails()
#' @name cumulative_tails



#' @rdname cumulative_tails
#' @export
calc_cumulative_tails <- function(object, ...) {
  UseMethod("calc_cumulative_tails")
}



#' @describeIn cumulative_tails Calculate the cumulative tails and saves to
#'   cevodata$models$cumulative_tails tibble
#' @export
calc_cumulative_tails.cevodata <- function(object,
                                           which_snvs = default_SNVs(object),
                                           column = get_frequency_measure_name(object, which_snvs),
                                           bins = 100,
                                           verbose = get_cevomod_verbosity(),
                                           ...) {
  cumulative_tails <- SNVs(object, which_snvs) |>
    calc_cumulative_tails(column = column, bins = bins, verbose = verbose)
  object$models[["cumulative_tails"]] <- cumulative_tails
  object
}



#' @describeIn cumulative_tails Calculate the cumulative tails
#' @export
calc_cumulative_tails.cevo_snvs <- function(object,
                                            column = get_frequency_measure_name(object),
                                            bins = 100,
                                            verbose = get_cevomod_verbosity(),
                                            ...) {
  msg("Calculating cumulative tails, using ", column, " column", verbose = verbose)

  res <- cut_f_intervals(object, column = column, bins = bins) |>
    filter(.data$f > 0) |>
    mutate(x = get_interval_centers(.data$f_interval), .after = "f_interval") |>
    group_by(.data$sample_id) |>
    .calc_cumulative_tails(bins) |>
    rename(f = "x", f_interval = "x_interval") |>
    ungroup()
  class(res) <- c("cevo_cumulative_tails_tbl", class(res))
  res
}



.calc_cumulative_tails <- function(dt, bins = 100) {
  breaks <- c(-1/bins, seq(0, 1, length.out = bins + 1))
  dt %>%
    mutate(
      x_interval = cut(.data$x, breaks = breaks),
      x_interval = as.character(.data$x_interval)
    ) %>%
    group_by(.data$x_interval, .add = TRUE) %>%
    summarise(n = n()) %>%
    mutate(x = get_interval_centers(.data$x_interval), .after = "x_interval") |>
    arrange(desc(.data$x), .by_group = TRUE) %>%
    mutate(
      y = cumsum(.data$n),
      y_scaled = round(.data$y / max(.data$y), digits = 4)
    )
}



#' @rdname cumulative_tails
#' @export
plot_cumulative_tails <- function(object, ...) {
  UseMethod("plot_cumulative_tails")
}



#' @describeIn cumulative_tails Shortcut to plot cum tails from SNVs dataframe
#'
#' @param scale_y scale y vaules to 1?
#' @param scales loglog/semilog
#' @param ... passed to geom_line()
#' @return ggplot obj
#' @export
plot_cumulative_tails.cevodata <- function(object, scale_y = TRUE, scales = "loglog", ...) {
  if (is.null(object$models$cumulative_tails)) {
    object <- calc_cumulative_tails(object)
  }
  object$models$cumulative_tails |>
    left_join(object$metadata, by = "sample_id") |>
    plot(scale_y = scale_y, scales = scales, ...)
}


#' Plot the cumulative tails
#'
#' @param x tibble with calc_cumulative_tails() results
#' @param scale_y logical
#' @param scales loglog/semilog
#' @param ... required by Generic
#' @return ggplot obj
#' @export
#'
#' @examples
#' data("test_data")
#' test_data |>
#'   calc_cumulative_tails() |>
#'   plot_cumulative_tails()
plot.cevo_cumulative_tails_tbl <- function(x, scale_y = TRUE, scales = "loglog", ...) {
  y <- if (scale_y) "y_scaled" else "y"

  plot_scales <- list(
    scale_y_log10(),
    scale_x_log10()
  )
  if (scales != "loglog") {
    plot_scales[[2]] <- NULL
  }

  x %>%
    ggplot(aes(.data$f, !!sym(y), color = .data$sample_id)) +
    geom_line(...) +
    plot_scales +
    theme_minimal() +
    labs(title = "Cumulative tails")
}


#' Plot Cumulative Tail
#'
#' Nice extention of the default ggplot2 stats. However, it should be easier to
#' use the plot_cumulative_tails() function.
#'
#' @inherit ggplot2::geom_histogram
#' @param bins number of bins
#' @param geom geom
#'
#' @examples
#' library(ggplot2)
#' data("tcga_brca_test")
#' snvs <- SNVs(tcga_brca_test)
#'
#' ggplot(snvs, aes(VAF, color = sample_id)) +
#'  stat_cumulative_tail()
#'
#' ggplot(snvs, aes(VAF, y = stat(y), color = sample_id)) +
#'  stat_cumulative_tail()
#' @export
stat_cumulative_tail <- function(mapping = NULL, data = NULL, geom = "point",
                                 position = "identity", na.rm = FALSE, show.legend = NA,
                                 inherit.aes = TRUE, bins = 100, ...) {
  layer(
    stat = StatCumulativeTail, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(bins = bins, na.rm = na.rm, ...)
  )
}


StatCumulativeTail <- ggproto("StatCumulativeTail", Stat,
  required_aes = c("x"),
  default_aes = aes(x = stat(x), y = stat(y_scaled)),
  compute_group = function(data, scales, bins = 100) {
    .calc_cumulative_tails(data, bins)
  }
)
