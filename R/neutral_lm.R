
#' Fitting neutral models
#'
#' Creates  cevodata$models$neutral_lm
#'
#' @param object SNVs tibble object
#' @param rsq_treshold R-squared tresholds to keep model as neutral
#' @param ... other arguments
#' @examples
#' data("tcga_brca_test")
#' snvs <- SNVs(tcga_brca_test) |>
#'   dplyr::filter(sample_id %in% c("TCGA-AC-A23H-01","TCGA-AN-A046-01"))
#'
#' cd <- init_cevodata("Test") |>
#'   add_SNV_data(snvs) |>
#'   calc_Mf_1f() |>
#'   calc_SFS() |>
#'   fit_neutral_lm(rsq_treshold = 0.99)
#'
#' plot(cd$models$Mf_1f, from = 0.05, to = 0.4, scale = FALSE) +
#'   layer_lm_fits(cd)
#' @name neutral_lm


#' @rdname neutral_lm
#' @export
fit_neutral_lm <- function(object, ...) {
  UseMethod("fit_neutral_lm")
}


#' @describeIn neutral_lm Fit Williams neutral models to the data
#' @export
fit_neutral_lm.cevodata <- function(object, rsq_treshold = 0.98, ...) {
  if (is.null(object$models$Mf_1f)) {
    stop("Run calc_Mf_1f() first!")
  }
  Mf_1f <- object$models$Mf_1f
  res <- Mf_1f |>
    filter(.data$`1/f` < Inf, .data$VAF < 0.5) |>
    select(.data$sample_id, .data$VAF, .data$`M(f)`, .data$`1/f`) |>
    nest(data = c(.data$VAF, .data$`M(f)`, .data$`1/f`)) |>
    mutate(fits = map(.data$data, fit_optimal_lm, rsq_treshold)) |>
    select(-.data$data) |>
    unnest(.data$fits)
  class(res) <- c("cevo_lm_models_tbl", class(res))
  object$models[["neutral_lm"]] <- res
  object <- calc_residuals(object)
  object
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


calc_residuals <- function(object, ...) {
  UseMethod("calc_residuals")
}

calc_residuals.cevodata <- function(object, ...) {
  neutral_lm <- object$models[["neutral_lm"]]
  sfs <- object$models[["SFS"]]
  if (is.null(neutral_lm) || is.null(sfs)) {
    stop("Calc SFS and and fit neutral lm first!")
  }

  binwidth <- get_average_interval(sfs$VAF)
  exp <- neutral_lm |>
    filter(.data$best) |>
    calc_powerlaw_curve(binwidth = binwidth) |>
    mutate(VAF = as.character(.data$f))
  sfs <- mutate(sfs, VAF = as.character(.data$VAF))

  dt <- exp |>
    select(.data$sample_id, .data$VAF, .data$neutral_pred) |>
    left_join(sfs, by = c("sample_id", "VAF"))

  residuals <- dt |>
    select(-.data$y_scaled) |>
    mutate(
      VAF = parse_double(.data$VAF),
      neutral_resid = .data$neutral_pred - .data$y,
      neutral_resid_clones = if_else(.data$neutral_resid > 0, 0, -.data$neutral_resid),
      sampling_rate = .data$neutral_resid / .data$neutral_pred
    )

  object$models[["residuals"]] <- residuals
  object
}


calc_powerlaw_curve <- function(lm_models, binwidth) {
  n_bins <- 1 / binwidth
  lm_models |>
    expand_grid(f = seq(0.01, 1, by = binwidth)) |>
    mutate(
      neutr = (.data$f >= .data$from & .data$f <= .data$to),
      neutral_pred = -(.data$a / n_bins) / .data$f^2
    )
}


get_average_interval <- function(vec) {
  vec |>
    unique() |>
    sort() |>
    diff() |>
    mean()
}


#' layer_neutral_tail
#'
#' @param object object
#' @param ... other arguments
#' @export
layer_neutral_tail <- function(object, ...) {
  UseMethod("layer_neutral_tail")
}

#' @describeIn neutral_lm Add M(f) ~ 1/f models layer to M(f) ~ 1/f plot
#'
#' @param cd cevodata
#' @param ... other params passed to geom_segment()
#' @export
layer_lm_fits <- function(cd, ...) {
  geom_segment(
    aes(
      x = 1/.data$from,
      xend = 1/.data$to,
      y = 1/.data$from * .data$a + .data$b,
      yend = 1/.data$to * .data$a + .data$b,
      color = .data$sample_id
    ),
    size = 1,
    data = cd$models$neutral_lm |>
      filter(, .data$best) |>
      left_join(cd$metadata, by = "sample_id"),
    show.legend = FALSE,
    ...
  )
}
