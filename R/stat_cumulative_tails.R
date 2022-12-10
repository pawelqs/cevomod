
#' Cumulative tails
#'
#' Creates  cevodata$models$cumulative_tails tbl with the groupping variables and:
#'   - n column with the number of mutations in the VAF interval
#'   - x and y columns describing the cumulative tails
#'   - y_scaled with y values scaled to the range 0-1
#'
#' @param object SNVs tibble object
#' @param digits resolution of the cumulative tails calculation
#' @param ... other arguments
#' @examples
#' data("tcga_brca_test")
#' tcga_brca_test |>
#'   calc_cumulative_tails()
#'
#' tcga_brca_test |>
#'   plot_cumulative_tails()
#' @name cumulative_tails


#' @rdname cumulative_tails
#' @export
calc_cumulative_tails <- function(object, ...) {
  UseMethod("calc_cumulative_tails")
}

#' @describeIn cumulative_tails Calculate the cumulative tails
#' @export
calc_cumulative_tails.cevodata <- function(object, digits = 2, ...) {
  cumulative_tails <- SNVs(object) |>
    group_by(.data$sample_id) |>
    calc_cumulative_tails(digits) %>%
    ungroup()
  class(cumulative_tails) <- c("cevo_cumulative_tails_tbl", class(cumulative_tails))
  object$models[["cumulative_tails"]] <- cumulative_tails
  object
}


#' @describeIn cumulative_tails Calculate the cumulative tails
#' @export
calc_cumulative_tails.tbl_df <- function(object, digits = 2, ...) {
  res <- object %>%
    rename(x = "VAF") %>%
    .calc_cumulative_tails(digits) %>%
    rename(VAF = "x")
  class(res) <- c("cevo_cumulative_tails_tbl", class(res))
  res
}


.calc_cumulative_tails <- function(dt, digits = 2) {
  dt %>%
    mutate(x = round(.data$x, digits = digits)) %>%
    group_by(.data$x, .add = TRUE) %>%
    summarise(n = n()) %>%
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
#' data("tcga_brca_test")
#' SNVs(tcga_brca_test) |>
#'   dplyr::group_by(sample_id) |>
#'   calc_cumulative_tails() |>
#'   plot()
plot.cevo_cumulative_tails_tbl <- function(x, scale_y = TRUE, scales = "loglog", ...) {

  y <- if (scale_y) "y_scaled" else "y"

  plot_scales <- list(
    scale_y_log10(),
    scale_x_log10()
  )
  if (scales != "loglog")
    plot_scales[[2]] <- NULL

  x %>%
    ggplot(aes(.data$VAF, !!sym(y), color = .data$sample_id)) +
    geom_line() +
    plot_scales +
    theme_minimal() +
    labs(title = "Cumulative tails")
}


#' @describeIn cumulative_tails Shortcut to plot cum tails from SNVs dataframe
#'
#' @param scale_y scale y vaules to 1?
#' @param ... passed to stat_cumulative_tail
#' @return ggplot obj
#' @export
plot_cumulative_tails.cevodata <- function(object, scale_y = TRUE, ...) {
  y <- NULL

  aes <- if (scale_y)
    aes(.data$VAF, color = .data$sample_id)
  else
    aes(.data$VAF, y = stat(y), color = .data$sample_id)

  SNVs(object) %>%
    left_join(object$metadata, by = "sample_id") |>
    ggplot(aes) +
    stat_cumulative_tail(...) +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal()
}


#' Plot Cumulative Tail
#'
#' @inherit ggplot2::geom_histogram
#' @param digits cumulative tail calculation avurracy
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
                                 inherit.aes = TRUE, digits = 2, ...) {
  layer(
    stat = StatCumulativeTail, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(digits = digits, na.rm = na.rm, ...)
  )
}


StatCumulativeTail <- ggproto("StatCumulativeTail", Stat,

  required_aes = c("x"),
  default_aes = aes(x = stat(x), y = stat(y_scaled)),

  compute_group = function(data, scales, digits = 2) {
    .calc_cumulative_tails(data, digits)
  }
)
