
#' Plot neutral exponential curve on SFS plot
#'
#' @param lm_models model fits from fit_neutral_lm
#' @param start VAF value to start plotting model fit
#' @param end VAF value to end plotting model fit=
#' @param show.legend show legend?
#' @param color line color
#' @param size line size
#' @param ... other arguments passed to geoms
#' @export
layer_exponential_fits <- function(lm_models,
                                   start = NA, end = 1, show.legend = FALSE,
                                   color = "black", size = 1, ...) {
  exp <- lm_models |>
    slice(1) |>
    calc_exp_curve() |>
    filter(f >= from - 0.02, f <= end)

  if (!is.na(start)) {
    exp <- filter(exp, f >= start)
  }

  list(
    geom_line(
      aes(f, n),
      data = exp,
      color = color, size = size, linetype = "dashed", show.legend = show.legend,
      ...
    ),
    geom_line(
      aes(f, n),
      data = exp |> filter(neutr),
      color = color, size = size, show.legend = show.legend,
      ...
    )
  )
}


calc_exp_curve <- function(lm_models) {
  lm_models |>
    expand_grid(f = seq(0.01, 1, by = 0.01)) |>
    mutate(
      neutr = (f >= from & f <= to),
      n = -(a / 90) / f^2
    )
}
