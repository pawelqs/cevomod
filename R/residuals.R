
get_residuals <- function(cd, model = cd$active_model) {
  cd$residuals[[model]]
}


#' Plot model residuals
#'
#' @param object cevodata object
#' @param mapping mapping elements to overwrite the default mapping
#' @param geom geom to use
#' @param fit_clones plot clonal fits?
#' @param ... other parameters
#' @name plot_residuals


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
    left_join(object$metadata, by = "sample_id")
  binomial_model_fitted <- !is.null(residuals[["binom_pred"]])
  default_mapping <- aes(.data$VAF, .data$neutral_resid_clones, color = .data$sample_id)
  final_mapping <- join_aes(default_mapping, mapping)
  clones_fit <- if (fit_clones && binomial_model_fitted) {
    fit_mapping <- aes(.data$VAF, .data$binom_pred, color = .data$sample_id, group = .data$sample_id)
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
