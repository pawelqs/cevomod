# --------------------------------- Fit ----------------------------------------

#' Fitting neutral models
#'
#' Fits a power-law component of the model:
#' \deqn{y(f) = \frac{\mu}{\beta n} \frac{1}{f^2}}
#' where \eqn{ \mu/\beta } is the mutation rate per effective cell division, and
#' n is the number of bins in the spectrum. The power-law exponent of this
#' model equals 2, as expected by
#' [Williams et al. (2016)](https://doi.org/10.1038/ng.3489) and
#' [Durrett (2013)](https://doi.org/10.1214/11-aap824).
#' This model is valid under the assumptions of exponential population growth,
#' constant mutation rate, and absence of selectively advantageous micro-clones
#' desctibed by
#' [Tung and Durrett (2021)](https://doi.org/10.1371/journal.pcbi.1008701).
#'
#' This model is fitted using the statistic
#' \deqn{M(f) = \frac{\mu}{\beta} \frac{1}{f}} described by
#' [Williams et al. (2016)](https://doi.org/10.1038/ng.3489).
#' Using the equations described in
#' [Williams et al. (2018)](https://doi.org/10.1038/s41588-018-0128-6)
#' evolutionary parameters of detected subclones can be calculated.
#'
#' @param object SNVs tibble object
#' @param rsq_treshold R-squared tresholds to keep model as neutral
#' @param lm_length length of the linear model fits
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
#'   intervalize_mutation_frequencies() |>
#'   calc_Mf_1f() |>
#'   calc_SFS() |>
#'   fit_powerlaw_tail_fixed(rsq_treshold = 0.99)
#'
#' plot(cd$models$Mf_1f, from = 0.05, to = 0.4, scale = FALSE) +
#'   layer_lm_fits(cd)
#' @name powerlaw_fixed_model
NULL


#' @rdname powerlaw_fixed_model
#' @export
fit_powerlaw_tail_fixed <- function(object, ...) {
  UseMethod("fit_powerlaw_tail_fixed")
}


#' @rdname powerlaw_fixed_model
#' @param pct_left drop pct of the lowerst frequency variants to improve fit
#' @param pct_right drop pct of the highest frequency variants to improve fit
#' @export
fit_powerlaw_tail_fixed.cevodata <- function(object,
                                             rsq_treshold = 0.98,
                                             lm_length = 0.05,
                                             name = "powerlaw_fixed",
                                             pct_left = 0.05, pct_right = 0.95,
                                             verbose = get_verbosity(),
                                             ...) {
  msg("Fitting williams neutral models...", verbose = verbose)

  Mf_1f <- get_Mf_1f(object)
  bounds <- SNVs(object) |>
    get_f_range(pct_left = pct_left, pct_right = pct_right)
  data <- Mf_1f |>
    left_join(bounds, by = "sample_id") |>
    filter(.data$f > .data$lower_bound, .data$f < .data$higher_bound) |>
    select("sample_id", "f", "M(f)", "1/f") |>
    nest_by(.data$sample_id)

  pb <- if (verbose) progress_bar$new(total = nrow(data)) else NULL
  coefs <- data |>
    reframe(
      model = "powerlaw_fixed",
      component = "Neutral tail",
      fit_optimal_lm(.data$data, rsq_treshold, lm_length = lm_length, pb)
    )
  residuals <- coefs |>
    filter(.data$best) |>
    calc_powerlaw_model_residuals(sfs = get_SFS(object))
  info <- list(f_column = attr(Mf_1f, "f_column"))

  models <- lst(coefs, residuals, info)
  class(models) <- c("cevo_powerlaw_models", "list")
  add_models(object, models, name = name)
}


fit_optimal_lm <- function(data, rsq_treshold = 0.98, lm_length = 0.05, pb = NULL) {
  min_val <- min(data$f)
  max_val <- max(data$f)
  grid <- expand_grid(
      from = seq(min_val, max_val, by = 0.01),
      to = seq(min_val, max_val, by = 0.01)
    ) |>
    mutate(length = .data$to - .data$from)

  desired_length <- if ((max_val - min_val) >= lm_length) lm_length else max(grid$length)
  grid <- grid |>
    filter(near(.data$length, desired_length))

  grid$data <- pmap(grid, prepare_Mf_1f_data, data = data)
  grid$fits <- map(grid$data, ~tidy_lm(.x$`1/f`, .x$`M(f)`))
  if (!is.null(pb)) pb$tick()
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
    filter(.data$f >= from, .data$f <= to)
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
    rsquared = if (stats::var(y) == 0) NA_real_ else stats::cor(y, x) ^ 2
  )
  res
}

# ----------------------------------- Get --------------------------------------

#' Get model coefficients
#' @param object cevodata object
#' @param model_name Model name
#' @param best_only Return only best models?
#' @export
get_model_coefficients <- function(object,
                                   model_name = active_models(object),
                                   best_only = TRUE) {
  coefs <- get_models(object, model_name)$coefs
  if (best_only) filter(coefs, .data$best) else coefs
}


#' Get model residuals
#' @param object cevodata object
#' @param model_name Model name
#' @export
get_model_residuals <- function(object, model_name = active_models(object)) {
  get_models(object, model_name)$residuals
}


# ----------------------------------- Plot -------------------------------------

#' Plot M(f) ~ 1/f fits
#' @param object cevodata object
#' @param ... other params
#' @export
plot_Mf_1f_fits <- function(object, model_name = "powerlaw_fixed", ...) {
  coefs <- get_model_coefficients(object, model_name) |>
    join_metadata(object)
  xmax <- max(coefs$to)
  xmax <- max(xmax, 0.25)
  plot_Mf_1f(object, scale = FALSE, to = xmax, ...) +
    layer_lm_fits(object, alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 90))
}


#' @describeIn powerlaw_fixed_model Add M(f) ~ 1/f models layer to M(f) ~ 1/f plot
#'
#' @param cd cevodata
#' @param model_name modelname
#' @param ... other params passed to geom_segment()
#' @export
layer_lm_fits <- function(cd, model_name = "powerlaw_fixed", ...) {
  geom_segment(
    aes(
      x = 1/.data$from,
      xend = 1/.data$to,
      y = 1/.data$from * .data$A + .data$b,
      yend = 1/.data$to * .data$A + .data$b
    ),
    linewidth = 1,
    data = get_model_coefficients(cd, model_name) |>
      join_metadata(cd),
    show.legend = FALSE,
    ...
  )
}


#' Plot 'a' coefficients for all fitted neutral models
#' @param object cevodata object
#' @param model_name modelname
#' @param ... other parameters passed to geom
#' @export
plot_neutral_A_coefficients <- function(object, ...) {
  UseMethod("plot_neutral_A_coefficients")
}


#' @export
plot_neutral_A_coefficients <- function(object, model_name = "powerlaw_fixed", ...) {
  get_model_coefficients(object, model_name, best_only = FALSE) |>
    join_metadata(object) |>
    ggplot() +
    aes(x = .data$from, xend = .data$to, y = .data$A, yend = .data$A, color = .data$best) +
    geom_segment(...) +
    facet_wrap(~.data$sample_id, scales = "free")
}
