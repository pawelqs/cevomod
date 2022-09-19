
#' Add M(f) ~ 1/f models layer to M(f) ~ 1/f plot
#'
#' @param data tibble with fits from fit_neutral_partial_models()
#' @param ... other params passed to geom_segment()
#' @export
#' @examples
#' Mf_1f <- snvs_test |>
#'   dplyr::filter(sample_id %in% c("TCGA-AC-A23H-01","TCGA-AN-A046-01")) |>
#'   dplyr::group_by(sample_id) |>
#'   calc_Mf_1f()
#' models <- fit_neutral_lm(Mf_1f, rsq_treshold = 0.99)
#' plot(Mf_1f, from = 0.05, to = 0.4, scale = FALSE) +
#'   layer_lm_fits(models)
layer_lm_fits <- function(data, ...) {
  geom_segment(
    aes(
      x = 1/.data$from,
      xend = 1/.data$to,
      y = 1/.data$from * .data$a + .data$b,
      yend = 1/.data$to * .data$a + .data$b,
      color = .data$sample_id
    ),
    size = 1,
    data = filter(data, .data$best),
    show.legend = FALSE,
    ...
  )
}


#' Fit Williams neutral models to the data
#'
#' @param dt tibble with results from calc_Mf_1f()
#' @param rsq_treshold R-squared tresholds to keep model as neutral
#' @return cevo_lm_models_tbl
#' @export
#' @inherit layer_lm_fits examples
fit_neutral_lm <- function(dt, rsq_treshold = 0.98) {
  res <- dt |>
    filter(.data$`1/f` < Inf, .data$VAF < 0.5) |>
    select(.data$sample_id, .data$VAF, .data$`M(f)`, .data$`1/f`) |>
    nest(data = c(.data$VAF, .data$`M(f)`, .data$`1/f`)) |>
    mutate(fits = map(.data$data, fit_optimal_lm, rsq_treshold)) |>
    select(-.data$data) |>
    unnest(.data$fits)
  class(res) <- c("cevo_lm_models_tbl", class(res))
  res
}


fit_optimal_lm <- function(dt, rsq_treshold = 0.98) {
  grid <- expand_grid(
      from = seq(0.05, 0.20, by = 0.01),
      to = seq(0.15, 0.3, by = 0.01)
    ) |>
    mutate(length = .data$to - .data$from) |>
    filter(.data$length > 0.05)
  grid$data <- pmap(grid, function(from, to, ...) filter(dt, .data$VAF >= from, .data$VAF <= to))
  grid$fits <- map(grid$data, function(dt) tidy_lm(dt$`1/f`, dt$`M(f)`))
  grid |>
    select(-.data$data) |>
    unnest(.data$fits) |>
    filter(.data$rsquared > rsq_treshold) |>
    arrange(desc(.data$length)) |>
    mutate(best = (row_number() == 1)) |>
    slice_head(n = 3)
}


tidy_lm <- function(x, y) {
  fit <- stats::lm(y ~ x)
  res <- tibble(
    a = fit$coefficients[["x"]],
    b = fit$coefficients[["(Intercept)"]],
    rsquared = stats::cor(y, x) ^ 2
  )
  res
}
