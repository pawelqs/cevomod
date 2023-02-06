
#' Fitting neutral models
#'
#' Creates  cevodata$models$neutral_model
#'
#' @param object SNVs tibble object
#' @param rsq_treshold R-squared tresholds to keep model as neutral
#' @param name name in the models' slot
#' @param verbose verbose?
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
#'   fit_neutral_models(rsq_treshold = 0.99)
#'
#' plot(cd$models$Mf_1f, from = 0.05, to = 0.4, scale = FALSE) +
#'   layer_lm_fits(cd)
#' @name neutral_model
NULL


#' @rdname neutral_model
#' @export
fit_neutral_models <- function(object, ...) {
  UseMethod("fit_neutral_models")
}


#' @describeIn neutral_model Fit Williams neutral models to the data
#' @export
fit_neutral_models.cevodata <- function(object,
                                        rsq_treshold = 0.98,
                                        name = "williams_neutral",
                                        verbose = TRUE, ...) {
  msg("Fitting williams neutral models...", verbose = verbose)
  Mf_1f <- get_Mf_1f(object)
  bounds <- get_VAF_range(SNVs(object))

  data <- Mf_1f |>
    left_join(bounds, by = "sample_id") |>
    filter(.data$VAF > .data$lower_bound, .data$VAF < .data$higher_bound) |>
    select("sample_id", "VAF", "M(f)", "1/f") |>
    nest(data = c("VAF", "M(f)", "1/f"))

  models <- data |>
    mutate(
      model = "williams",
      component = "Neutral tail",
      fits = map(.data$data, fit_optimal_lm, rsq_treshold)
    ) |>
    select(-"data") |>
    unnest("fits")
  class(models) <- c("cevo_powerlaw_models", class(models))

  object$models[[name]] <- models
  object$active_models <- name
  object <- calc_powerlaw_model_residuals(object, models_name = name)
  object
}


fit_optimal_lm <- function(data, rsq_treshold = 0.98) {
  min_val <- min(data$VAF)
  max_val <- max(data$VAF)
  grid <- expand_grid(
      from = seq(min_val, max_val, by = 0.01),
      to = seq(min_val, max_val, by = 0.01)
    ) |>
    mutate(length = .data$to - .data$from) |>
    filter(near(.data$length, 0.05))
  grid$data <- pmap(grid, prepare_Mf_1f_data, data = data)
  grid$fits <- map(grid$data, ~tidy_lm(.x$`1/f`, .x$`M(f)`))
  grid |>
    select(-"data") |>
    unnest("fits") |>
    filter(.data$rsquared > rsq_treshold) |>
    rename(A = "a") |>
    mutate(alpha = 2, .before = "rsquared") |>
    arrange(.data$A) |>
    mutate(best = (row_number() == 1))
}


prepare_Mf_1f_data <- function(from, to, data, ...) {
  data |>
    filter(.data$VAF >= from, .data$VAF <= to)
    # Following may be used to fit the model exactly as Williams, without the
    # Intercept term (lm(y ~ x + 0)). It does not affect slope coefficient,
    # thus I keep the current implementation for visualization purposes
    #
    # mutate(
    #   y = `M(f)` - min(`M(f)`),
    #   x = `1/f` - min(`1/f`)
    # )
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


#' @rdname neutral_model
#' @export
get_neutral_models <- function(object, ...) {
  UseMethod("get_neutral_models")
}


#' @describeIn neutral_model Get neutral models
#' @param best_only return only the best fits
#' @export
get_neutral_models <- function(object, best_only = TRUE, ...) {
  models <- object$models$neutral_models
  if (best_only) {
    filter(models, .data$best)
  } else {
    models
  }
}


#' Plot M(f) ~ 1/f fits
#' @param object cevodata object
#' @param ... other params
#' @export
plot_Mf_1f_fits <- function(object, ...) {
  plot_Mf_1f(object, scale = FALSE, ...) +
    layer_lm_fits(object, alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 90))
}


#' @describeIn neutral_model Add M(f) ~ 1/f models layer to M(f) ~ 1/f plot
#'
#' @param cd cevodata
#' @param ... other params passed to geom_segment()
#' @export
layer_lm_fits <- function(cd, ...) {
  geom_segment(
    aes(
      x = 1/.data$from,
      xend = 1/.data$to,
      y = 1/.data$from * .data$A + .data$b,
      yend = 1/.data$to * .data$A + .data$b
    ),
    size = 1,
    data = get_neutral_models(cd) |>
      left_join(cd$metadata, by = "sample_id"),
    show.legend = FALSE,
    ...
  )
}


#' Plot 'a' coefficients for all fitted neutral models
#' @param object cevodata object
#' @param ... other parameters passed to geom
#' @export
plot_neutral_A_coefficients <- function(object, ...) {
  UseMethod("plot_neutral_A_coefficients")
}


#' @export
plot_neutral_A_coefficients <- function(object, ...) {
  get_neutral_models(object, best_only = FALSE) |>
    ggplot() +
    aes(x = .data$from, xend = .data$to, y = .data$A, yend = .data$A, color = .data$best) +
    geom_segment(...) +
    facet_wrap(~.data$sample_id, scales = "free") +
    theme_minimal() +
    scale_color_brewer(palette = "Dark2")
}
