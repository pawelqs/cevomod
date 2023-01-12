
#' Plot cevodata models
#' @param object cevodata object
#' @param neutral_tail TRUE,
#' @param binomial_layer FALSE,
#' @param subclones TRUE,
#' @param final_fit TRUE,
#' @param ... other arguments
#' @name plot_models


#' @rdname plot_models
#' @export
plot_models <- function(object, ...) {
  UseMethod("plot_models")
}


#' @rdname plot_models
#' @export
plot_models.cevodata <- function(object,
                                 neutral_tail = TRUE,
                                 binomial_layer = FALSE,
                                 subclones = TRUE,
                                 final_fit = TRUE,
                                 ...) {

  neutral_lm_fitted <- !is.null(object$models[["neutral_models"]])
  subclones_fitted <- !is.null(object$models[["binomial_models"]])

  neutral_models <- get_neutral_models(object) |>
    select("sample_id", "from", "to")

  resid <- get_residuals(object) |>
    left_join(neutral_models, by = "sample_id") |>
    group_by(.data$sample_id) |>
    mutate(
      ylim = max(.data$SFS) * 1.2,
      neutral_pred = if_else(.data$neutral_pred > .data$ylim, .data$ylim, .data$neutral_pred)
    ) |>
    ungroup() |>
    mutate(neutr = .data$VAF >= .data$from & .data$VAF <= .data$to)

  model_layers <- list(
    if (neutral_tail && neutral_lm_fitted) {
      geom_area(
        aes(.data$VAF, .data$neutral_pred),
        data = resid, # |> filter(.data$VAF >= 0),
        fill = "white", color = "gray90",
        alpha = 0.3,
        size = 0.5, show.legend = FALSE,
        stat = "identity"
      )
    },
    # if (neutral_tail && neutral_lm_fitted) {
    #   geom_line(
    #     aes(.data$VAF, .data$neutral_pred),
    #     data = resid |> filter(.data$neutr, .data$neutral_pred < .data$ylim),
    #     size = 1, show.legend = FALSE
    #   )
    # },
    if (binomial_layer && subclones_fitted) {
      geom_line(
        aes(.data$VAF, .data$binom_pred),
        data = resid,
        size = 1, color = "black", linetype = "dashed"
      )
    },
    if (subclones && subclones_fitted) {
      dt <- resid |>
        pivot_longer(
          cols = c("Clone", starts_with("Subclone")),
          names_to = "component",
          values_to = "pred"
        ) |>
        filter(!is.na(.data$pred))
      geom_area(
        aes(.data$VAF, .data$pred, group = .data$component),
        data = dt,
        position = "identity",
        size = 1, alpha = 0.3, color = "black", show.legend = FALSE
      )
    },
    if (final_fit && neutral_lm_fitted && subclones_fitted) {
      geom_line(
        aes(.data$VAF, .data$model_pred),
        data = resid |> filter(.data$neutral_pred < .data$ylim),
        size = 1, color = "red"
      )
    }
  )

  plot_SFS(object, geom = "bar") +
    model_layers +
    facet_wrap(~.data$sample_id, scales = "free_y")
}


plot_clones <- function(clones) {
  clones |>
    get_binomial_predictions() |>
    select(-"binom_pred") |>
    pivot_longer(-.data$VAF) |>
    ggplot(aes(.data$VAF, .data$value, group = .data$name)) +
    geom_point()
}

