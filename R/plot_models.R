
#' Plot cevodata models
#' @param object cevodata object
#' @param models_name models name
#' @param show_neutral_tail TRUE,
#' @param show_binomial_layer FALSE,
#' @param show_subclones TRUE,
#' @param show_final_fit TRUE,
#' @param ... other arguments passed to plot_SFS()
#' @name plot_models


#' @rdname plot_models
#' @export
plot_models <- function(object, ...) {
  UseMethod("plot_models")
}


#' @rdname plot_models
#' @param neutral_tail_alpha 0.3
#' @param neutral_tail_size 0.5
#' @param neutral_tail_fill "white"
#' @param neutral_tail_color "gray90"
#' @param binomial_layer_color "black"
#' @param final_fit_color "red"
#' @param final_fit_size 1
#' @param nrow passed to facet_wrap
#' @param ncol passed to facet_wrap
#' @export
plot_models.cevodata <- function(object,
                                 models_name = active_models(object),
                                 show_neutral_tail = TRUE,
                                 show_binomial_layer = FALSE,
                                 show_subclones = TRUE,
                                 show_final_fit = TRUE,
                                 neutral_tail_alpha = 0.3,
                                 neutral_tail_size = 0.5,
                                 neutral_tail_fill = "white",
                                 neutral_tail_color = "gray90",
                                 binomial_layer_color = "black",
                                 final_fit_color = "red",
                                 final_fit_size = 1,
                                 nrow = NULL, ncol = NULL,
                                 ...) {

  models <- get_models(object, models_name)
  neutral_lm_fitted <- "alpha" %in% names(models)
  subclones_fitted <- "cellularity" %in% names(models)

  resid <- get_residuals(object, models_name) |>
    left_join(object$metadata, by = "sample_id") |>
    mutate(sample_id = parse_factor(.data$sample_id, levels = object$metadata$sample_id)) |>
    group_by(.data$sample_id) |>
    mutate(
      ylim = max(.data$SFS) * 1.2,
      powerlaw_pred = case_when(
        .data$powerlaw_pred > .data$ylim ~ Inf,
        .data$VAF < 0 ~ Inf,
        TRUE ~ .data$powerlaw_pred
      )
    ) |>
    ungroup()

  model_layers <- list(
    if (show_neutral_tail && neutral_lm_fitted) {
      geom_area(
        aes(.data$VAF, .data$powerlaw_pred),
        data = resid,
        fill = neutral_tail_fill, color = neutral_tail_color,
        alpha = neutral_tail_alpha,
        size = neutral_tail_size, show.legend = FALSE,
        stat = "identity"
      )
    },
    if (show_binomial_layer && subclones_fitted) {
      geom_line(
        aes(.data$VAF, .data$binom_pred),
        data = resid,
        size = 1, color = binomial_layer_color, linetype = "dashed"
      )
    },
    if (show_subclones && subclones_fitted) {
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
    if (show_final_fit && neutral_lm_fitted && subclones_fitted) {
      geom_line(
        aes(.data$VAF, .data$model_pred),
        data = resid |> filter(.data$powerlaw_pred < .data$ylim),
        size = final_fit_size, color = final_fit_color
      )
    }
  )

  plot_SFS(object, geom = "bar", ...) +
    model_layers +
    facet_wrap(~.data$sample_id, scales = "free_y", nrow = nrow, ncol = ncol)
}


plot_clones <- function(clones) {
  clones |>
    get_binomial_predictions() |>
    select(-"binom_pred") |>
    pivot_longer(-.data$VAF) |>
    ggplot(aes(.data$VAF, .data$value, group = .data$name)) +
    geom_point()
}


#' Plot power-law curve
#' @param A A
#' @param alpha power
#' @param mapping mapping, x is required
#' @param ylim max y allowed
#' @param color color
#' @param ... other arguments passed to geom_line
#' @export
geom_powerlaw <- function(A, alpha, mapping, ylim = 1000, color = "#54b4FA", ...) {
  . <- NULL
  geom_line(
    join_aes(aes(.data$VAF, .data$y), mapping),
    data = . %>%
      mutate(y = A * 1 / .data$VAF ^ alpha) %>%
      filter(.data$y <= ylim, .data$y > 0),
    linewidth = 1.5,
    color = color,
    ...
  )
}


#' Compare fits from many models
#' @param object cevodata object
#' @param model_names models to compare
#' @param column_name residuals_* column to plot
#' @param linetype solid
#' @param linewidth 1
#' @param ... other arguments passed to plot_SFS()
#' @export
compare_models <- function(object, model_names, column_name,
                           linetype = "solid", linewidth = 1, ...) {
  ylimits <- get_SFS(object) |>
    group_by(.data$sample_id) |>
    summarise(ylim = 1.5 * max(.data$y))

  resids <- model_names %>%
    set_names(model_names) |>
    map(~get_residuals(object, .x)) |>
    bind_rows(.id = "model_name") |>
    left_join(ylimits, by = "sample_id") |>
    left_join(object$metadata, by = "sample_id") |>
    filter(!!sym(column_name) < .data$ylim, .data$VAF >= 0)

  plot_SFS(object, geom = "bar", ...) +
    geom_line(
      aes(y = !!sym(column_name), color = .data$model_name, group = .data$model_name),
      data = resids,
      linetype = linetype, linewidth = linewidth
    ) +
    facet_wrap(~.data$sample_id, scales = "free_y")
}

