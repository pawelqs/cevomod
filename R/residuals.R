
#' Get model residuals
#' @param cd cevodata object
#' @param models_name name of the models
#' @export
get_residuals <- function(cd, models_name = cd$active_model) {
  slot_name <- paste0("residuals_", models_name)
  residuals <- cd$misc[[slot_name]]
  if (is.null(residuals)) {
    stop(slot_name, "slot empty. Fit apropriate model first!")
  }
  residuals
}


calc_powerlaw_model_residuals <- function(powerlaw_coefs, sfs, ...) {
  optional_cols <- c("from", "to", "b") |> intersect(colnames(powerlaw_coefs))
  from_to_cols_present <- all(c("from", "to") %in% optional_cols)
  powerlaw_coefs <- powerlaw_coefs |>
    filter(!is.na(.data$A), !is.na(.data$alpha)) |>
    select("sample_id", any_of("resample_id"), "A", "alpha", all_of(optional_cols))
  nbins <- summarise(sfs, nbins = n() - 1, .by = "sample_id") # zero bin does not count

  residuals <- sfs |>
    select("sample_id", "f_interval", "f", SFS = "y") |>
    inner_join(powerlaw_coefs, by = "sample_id") |>
    left_join(nbins, by = "sample_id") |>
    mutate(
      neutr = if (from_to_cols_present) {
        .data$f >= .data$from & .data$f <= .data$to
      } else NA,
    ) |>
    mutate(
      powerlaw_pred = calc_powerlaw_curve(.data$f, .data$A, .data$alpha, .data$nbins),
      powerlaw_resid = .data$powerlaw_pred - .data$SFS,
      powerlaw_resid_clones = if_else(.data$powerlaw_resid > 0, 0, -.data$powerlaw_resid),
      sampling_rate = .data$powerlaw_resid / .data$powerlaw_pred,
      model_resid = .data$powerlaw_resid,
    ) |>
    select(-("nbins":"alpha"))
  class(residuals) <- class(tibble())
  residuals
}


calc_powerlaw_curve <- function(f, A, alpha, nbins) {
  if_else(f < 0, 0, (A / nbins) / f^alpha)
}


#' Plot model residuals
#'
#' @param object cevodata object
#' @param mapping mapping elements to overwrite the default mapping
#' @param geom geom to use
#' @param models_name models_name
#' @param fit_clones plot clonal fits?
#' @param ... other parameters
#' @name plot_residuals
NULL



#' @describeIn plot_residuals Plot sampling rate
#' @export
plot_sampling_rate <- function(object, mapping = NULL, geom = geom_point, ...) {
  residuals <- get_residuals(object) |>
    left_join(object$metadata, by = "sample_id") |>
    filter(.data$f >= 0)
  default_mapping <- aes(.data$f, .data$sampling_rate, color = .data$sample_id)
  final_mapping <- join_aes(default_mapping, mapping)
  ggplot(residuals) +
    geom(...) +
    final_mapping +
    coord_cartesian(xlim = c(0, .2), ylim = c(0, 1)) +
    labs(y = "Sampling Rate") +
    theme_minimal()
}



#' @describeIn plot_residuals Plot residuals of the neutral model
#' @export
plot_residuals_powerlaw_model <- function(object,
                                          models_name = active_models(object),
                                          mapping = NULL,
                                          geom = geom_point,
                                          fit_clones = TRUE,
                                          ...) {
  residuals <- get_residuals(object, models_name) |>
    left_join(object$metadata, by = "sample_id") |>
    group_by(.data$sample_id) |>
    mutate(width = 0.9 / n())
  binomial_model_fitted <- !is.null(residuals[["binom_pred"]])
  default_mapping <- aes(.data$f, .data$powerlaw_resid_clones, group = .data$sample_id, width = .data$width)
  final_mapping <- join_aes(default_mapping, mapping)
  clones_fit <- if (fit_clones && binomial_model_fitted) {
    fit_mapping <- aes(.data$f, .data$binom_pred, group = .data$sample_id)
    final_fit_mapping <- join_aes(fit_mapping, mapping)
    geom_line(final_fit_mapping, color = "black")
  }
  ggplot(residuals) +
    geom(...) +
    clones_fit +
    final_mapping +
    labs(y = "Residuals") +
    theme_minimal()
}


#' @describeIn plot_residuals Plot residuals of the full model
#' @export
plot_residuals_full_model <- function(object,
                                      mapping = NULL,
                                      geom = geom_point,
                                      ...) {
  residuals <- get_residuals(object) |>
    left_join(object$metadata, by = "sample_id")
  default_mapping <- aes(.data$f, .data$model_resid, color = .data$sample_id)
  final_mapping <- join_aes(default_mapping, mapping)
  ggplot(residuals) +
    geom(...) +
    final_mapping +
    labs(y = "Residuals") +
    coord_cartesian(ylim = c(-1000, NA_integer_)) +
    theme_minimal()
}


#' @describeIn plot_residuals Plot binomial fits vs powerlaw residuals (barplot)
#' @export
plot_binomial_fits_vs_powerlaw_residuals_bars <- function(
        object,
        models_name = active_models(object),
        mapping = NULL,
        geom = geom_bar,
        fit_clones = TRUE,
        ...) {
  residuals <- get_residuals(object, models_name)
  binomial_model_fitted <- !is.null(residuals[["binom_pred"]])
  if (!binomial_model_fitted) {
    stop("Fit subclones first!")
  }

  dt <- residuals |>
    select("sample_id", "f", resid = "powerlaw_resid_clones", pred = "binom_pred") |>
    group_by(.data$sample_id) |>
    mutate(width = 0.9 / n()) |>
    pivot_longer(
      -c("sample_id", "f", "width"),
      names_to = "variable", values_to = "value"
    ) |>
    left_join(object$metadata, by = "sample_id")

  ggplot(dt) +
    aes(.data$f, .data$value, width = .data$width, fill = .data$variable) +
    geom_bar(stat = "identity", position = "identity", alpha = 0.5) +
    facet_wrap(~.data$sample_id, scales = "free")
}
