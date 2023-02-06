
get_residuals <- function(cd, model = cd$active_model) {
  slot_name <- paste0("residuals_", model)
  residuals <- cd$misc[[slot_name]]
  if (is.null(residuals)) {
    stop(slot_name, "slot empty. Fit apropriate model first!")
  }
  residuals
}


calc_powerlaw_model_residuals <- function(object, models_name, ...) {
  powerlaw_models <- get_models(object, models_name)
  optional_cols <- c("from", "to", "b") |> intersect(colnames(powerlaw_models))
  from_to_cols_present <- all(c("from", "to") %in% optional_cols)
  powerlaw_models <- powerlaw_models |>
    select("sample_id", "A", "alpha", all_of(optional_cols))
  sfs <- get_SFS(object)
  nbins <- get_sample_sequencing_depths(SNVs(object)) |>
    transmute(.data$sample_id, nbins = .data$median_DP)

  residuals <- sfs |>
    select("sample_id", "VAF_interval", "VAF", SFS = "y") |>
    inner_join(powerlaw_models, by = "sample_id") |>
    left_join(nbins, by = "sample_id") |>
    mutate(
      neutr = if (from_to_cols_present) {
        .data$VAF >= .data$from & .data$VAF <= .data$to
        } else NA,
    ) |>
    mutate(
      neutral_pred = calc_powerlaw_curve(.data$VAF, .data$A, .data$alpha, .data$nbins),
      neutral_resid = .data$neutral_pred - .data$SFS,
      neutral_resid_clones = if_else(.data$neutral_resid > 0, 0, -.data$neutral_resid),
      sampling_rate = .data$neutral_resid / .data$neutral_pred,
      model_resid = .data$neutral_resid,
    ) |>
    select(-("nbins":"alpha"))

  slot_name <- paste0("residuals_", models_name)
  object$misc[[slot_name]] <- residuals
  object
}


calc_powerlaw_curve <- function(VAF, A, alpha, nbins) {
  if_else(VAF < 0, 0, (A / nbins) / VAF^alpha)
}


#' Plot model residuals
#'
#' @param object cevodata object
#' @param mapping mapping elements to overwrite the default mapping
#' @param geom geom to use
#' @param fit_clones plot clonal fits?
#' @param ... other parameters
#' @name plot_residuals
NULL


#' @rdname plot_residuals
#' @export
plot_sampling_rate <- function(object, ...) {
  UseMethod("plot_sampling_rate")
}


#' @describeIn plot_residuals Plot sampling rate
#' @export
plot_sampling_rate.cevodata <- function(object, mapping = NULL, geom = geom_point, ...) {
  residuals <- get_residuals(object) |>
    left_join(object$metadata, by = "sample_id") |>
    filter(.data$VAF >= 0)
  default_mapping <- aes(.data$VAF, .data$sampling_rate, color = .data$sample_id)
  final_mapping <- join_aes(default_mapping, mapping)
  ggplot(residuals) +
    geom(...) +
    final_mapping +
    coord_cartesian(xlim = c(0, .2), ylim = c(0, 1)) +
    labs(y = "Sampling Rate") +
    theme_minimal()
}


#' @rdname plot_residuals
#' @export
plot_residuals_neutral_model <- function(object, ...) {
  UseMethod("plot_residuals_neutral_model")
}


#' @describeIn plot_residuals Plot residuals of the neutral model
#' @export
plot_residuals_neutral_model.cevodata <- function(object,
                                                  mapping = NULL,
                                                  geom = geom_point,
                                                  fit_clones = TRUE,
                                                  ...) {
  residuals <- get_residuals(object) |>
    left_join(object$metadata, by = "sample_id") |>
    group_by(.data$sample_id) |>
    mutate(width = 0.9 / n())
  binomial_model_fitted <- !is.null(residuals[["binom_pred"]])
  default_mapping <- aes(.data$VAF, .data$neutral_resid_clones, group = .data$sample_id, width = .data$width)
  final_mapping <- join_aes(default_mapping, mapping)
  clones_fit <- if (fit_clones && binomial_model_fitted) {
    fit_mapping <- aes(.data$VAF, .data$binom_pred, group = .data$sample_id)
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


#' @rdname plot_residuals
#' @export
plot_residuals_full_model <- function(object, ...) {
  UseMethod("plot_residuals_full_model")
}


#' @describeIn plot_residuals Plot residuals of the full model
#' @export
plot_residuals_full_model.cevodata <- function(object,
                                               mapping = NULL,
                                               geom = geom_point,
                                               ...) {
  residuals <- get_residuals(object) |>
    left_join(object$metadata, by = "sample_id")
  default_mapping <- aes(.data$VAF, .data$model_resid, color = .data$sample_id)
  final_mapping <- join_aes(default_mapping, mapping)
  # y_min <- residuals |>
  #   filter(VAF > .25) |>
  #   pull(model_resid) |>
  #   min()
  ggplot(residuals) +
    geom(...) +
    final_mapping +
    labs(y = "Residuals") +
    coord_cartesian(ylim = c(-1000, NA_integer_)) +
    theme_minimal()
}
