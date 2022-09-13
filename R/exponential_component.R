
#' Plot neutral exponential curve on SFS plot
#'
#' @param models model fits from fit_neutral_lm
#' @param start VAF value to start plotting model fit
#' @param end VAF value to end plotting model fit
#' @param color line color
#' @param size line size
#' @param ... other arguments passed to geoms
#' @export
layer_exponential_fits <- function(models,
                                   start = 0.08, end = 1, show.legend = FALSE,
                                   color = "black", size = 1, ...) {
  exp <- models |>
    slice(1) |>
    expand_grid(x = seq(start, end, by = 0.01)) |>
    mutate(
      neutr = (x >= from & x <= to),
      y = -(a / 90) / x^2
    )

  list(
    geom_line(
      aes(x, y),
      data = exp,
      color = color, size = size, linetype = "dashed", show.legend = show.legend,
      ...
    ),
    geom_line(
      aes(x, y),
      data = exp |> filter(neutr),
      color = color, size = size, show.legend = show.legend,
      ...
    )
  )
}
