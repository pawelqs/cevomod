
#' @describeIn layer_neutral_tail Plot neutral exponential curve on SFS plot
#'
#' @param object model fits from fit_neutral_lm
#' @param start VAF value to start plotting model fit
#' @param end VAF value to end plotting model fit
#' @param binwidth binwidth, necessary to scale the model parameter so it can
#'   fit the binarized SFS correctly
#' @param show.legend show legend?
#' @param color line color
#' @param size line size
#' @param ... other arguments passed to geoms
#' @export
layer_neutral_tail.cevodata <- function(object,
                                        start = NA, end = 1,
                                        binwidth = 0.01,
                                        show.legend = FALSE,
                                        color = "black", size = 1, ...) {
  exp <- object$models$neutral_lm |>
    filter(.data$best) |>
    calc_powerlaw_curve(binwidth) |>
    filter(.data$f >= .data$from - 0.02, .data$f <= end) |>
    left_join(object$metadata, by = "sample_id")

  if (!is.na(start)) {
    exp <- filter(exp, .data$f >= start)
  }

  list(
    geom_line(
      aes(.data$f, .data$neutral_pred),
      data = exp,
      color = color, size = size, linetype = "dashed", show.legend = show.legend,
      ...
    ),
    geom_line(
      aes(.data$f, .data$neutral_pred),
      data = exp |> filter(.data$neutr),
      color = color, size = size, show.legend = show.legend,
      ...
    )
  )
}


