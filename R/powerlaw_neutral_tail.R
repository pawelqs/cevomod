
#' Plot neutral exponential curve on SFS plot
#'
#' @param lm_models model fits from fit_neutral_lm
#' @param start VAF value to start plotting model fit
#' @param end VAF value to end plotting model fit
#' @param binwidth binwidth, necessary to scale the model parameter so it can
#'   fit the binarized SFS correctly
#' @param show.legend show legend?
#' @param color line color
#' @param size line size
#' @param ... other arguments passed to geoms
#' @export
layer_neutral_tail <- function(lm_models,
                               start = NA, end = 1,
                               binwidth = 0.01,
                               show.legend = FALSE,
                               color = "black", size = 1, ...) {
  exp <- lm_models |>
    slice(1) |>
    calc_powerlaw_curve(binwidth) |>
    filter(.data$f >= .data$from - 0.02, .data$f <= end)

  if (!is.na(start)) {
    exp <- filter(exp, .data$f >= start)
  }

  list(
    geom_line(
      aes(.data$f, .data$n),
      data = exp,
      color = color, size = size, linetype = "dashed", show.legend = show.legend,
      ...
    ),
    geom_line(
      aes(.data$f, .data$n),
      data = exp |> filter(.data$neutr),
      color = color, size = size, show.legend = show.legend,
      ...
    )
  )
}


calc_powerlaw_curve <- function(lm_models, binwidth) {
  n_bins <- 1 / binwidth
  lm_models |>
    expand_grid(f = seq(0.01, 1, by = binwidth)) |>
    mutate(
      neutr = (.data$f >= .data$from & .data$f <= .data$to),
      n = -(.data$a / n_bins) / .data$f^2
    )
}


#' Estimate sampling rate
#'
#' Uses experimental SFS and power-law model to estimate the sampling rate
#'
#' @param sfs SFS
#' @param lm_models output from fit_neutral_lm()
#'
#' @return tibble
#' @export
estimate_sampling_rate <- function(sfs, lm_models) {
  binwidth <- get_average_interval(sfs$VAF)
  exp <- calc_powerlaw_curve(lm_models, binwidth = binwidth) |>
    mutate(VAF = as.character(.data$f))
  sfs <- mutate(sfs, VAF = as.character(VAF))

  patient_id_present <- "patient_id" %in% names(sfs)
  dt <- if (patient_id_present) {
    exp |>
      select(.data$patient_id, .data$sample_id, .data$VAF, .data$n) |>
      left_join(sfs, by = c("patient_id", "sample_id", "VAF"))
  } else {
    exp |>
      select(.data$sample_id, .data$VAF, .data$n) |>
      left_join(sfs, by = c("sample_id", "VAF"))
  }

  sampling_stats <- dt |>
    select(-.data$y_scaled) |>
    filter(.data$VAF < 0.2) |>
    mutate(
      VAF = parse_double(.data$VAF),
      err = .data$n - .data$y,
      sampling_rate = .data$err / .data$n
    )
  sampling_stats
}


get_average_interval <- function(vec) {
  vec |>
    unique() |>
    sort() |>
    diff() |>
    mean()
}
