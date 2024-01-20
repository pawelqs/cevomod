
## ------------------------- plot_models() function ---------------------------

#' Plot cevodata models
#' @param object `<cevodata>` object
#' @param models_name Name of the models to plot
#' @param show_neutral_tail `<lgl>`
#' @param show_binomial_layer `<lgl>`
#' @param show_subclones `<lgl>`
#' @param show_final_fit `<lgl>`
#' @param ... Other arguments passed to plot_SFS()
#' @name plot_models



#' @rdname plot_models
#' @export
plot_models <- function(object, ...) {
  UseMethod("plot_models")
}



#' @rdname plot_models
#' @param params_neutral_tail List of non-default params for neutral tail geom
#' @param params_bootstraps List of non-default params for neutral tail bootstraps
#'   geom
#' @param params_binomial List of non-default params for binomial geom
#' @param params_subclones List of non-default params for subclones geom
#' @param params_final_fit List of non-default params for final fit geom
#' @param nrow Passed to facet_wrap
#' @param ncol Passed to facet_wrap
#' @export
plot_models.cevodata <- function(object,
                                 models_name = active_models(object),
                                 show_neutral_tail = TRUE,
                                 show_binomial_layer = FALSE,
                                 show_subclones = TRUE,
                                 show_final_fit = TRUE,
                                 params_neutral_tail = list(),
                                 params_bootstraps = list(),
                                 params_binomial = list(),
                                 params_subclones = list(),
                                 params_final_fit = list(),
                                 nrow = NULL, ncol = NULL,
                                 ...) {
  models <- get_models(object, models_name)
  neutral_lm_fitted <- "alpha" %in% names(models)
  subclones_fitted <- "cellularity" %in% names(models)
  bootstraped <- "cevo_bootstrap_powerlaw_models" %in% class(models)

  resid <- get_residuals(object, models_name) |>
    left_join(object$metadata, by = "sample_id") |>
    mutate(sample_id = parse_factor(.data$sample_id, levels = object$metadata$sample_id)) |>
    trim_powerlaw_pred()

  model_layers <- list(
    if (show_neutral_tail && neutral_lm_fitted && !bootstraped) {
      geom_area_neutral_tail(resid, params_neutral_tail)
    },
    if (show_neutral_tail && neutral_lm_fitted && bootstraped) {
      geom_line_bootstraps(resid, params_bootstraps)
    },
    if (show_binomial_layer && subclones_fitted) {
      geom_line_binomial(resid, params_binomial)
    },
    if (show_subclones && subclones_fitted) {
      geom_are_subclones(resid, params_subclones)
    },
    if (show_final_fit && neutral_lm_fitted && subclones_fitted) {
      geom_line_final_fit(resid, params_final_fit)
    }
  )

  plot_SFS(object, geom = "bar", ...) +
    model_layers +
    facet_wrap(~.data$sample_id, scales = "free_y", nrow = nrow, ncol = ncol)
}



trim_powerlaw_pred <- function(resid, limit = 1.2) {
  resid |>
    group_by(.data$sample_id) |>
    mutate(
      ylim = max(.data$SFS) * limit,
      powerlaw_pred = case_when(
        .data$powerlaw_pred > .data$ylim ~ Inf,
        .data$f < 0 ~ Inf,
        TRUE ~ .data$powerlaw_pred
      )
    ) |>
    ungroup()
}


geom_area_neutral_tail <- function(resid, params = list()) {
  geom_area(
    aes(.data$f, .data$powerlaw_pred),
    data = resid,
    fill = if (is.null(params$fill)) "white" else params$fill,
    color = if (is.null(params$color)) "gray90" else params$color,
    alpha = if (is.null(params$alpha)) 0.3 else params$alpha,
    size = if (is.null(params$size)) 0.5 else params$size,
    show.legend = FALSE,
    stat = "identity"
  )
}



geom_line_bootstraps <- function(resid, params = list()) {
  geom_line(
    aes(.data$f, .data$powerlaw_pred, group = .data$resample_id),
    data = resid,
    alpha = if (is.null(params$alpha)) 1 else params$alpha,
    color = if (is.null(params$color)) "white" else params$color,
    linetype = if (is.null(params$linetype)) "solid" else params$linetype,
    linewidth = if (is.null(params$linewidth)) 0.3 else params$linewidth,
    show.legend = FALSE,
    stat = "identity"
  )
}



geom_line_binomial <- function(resid, params = list()) {
  geom_line(
    aes(.data$f, .data$binom_pred),
    data = resid,
    alpha = if (is.null(params$alpha)) 1 else params$alpha,
    color = if (is.null(params$color)) "black" else params$color,
    linetype = if (is.null(params$linetype)) "dashed" else params$linetype,
    linewidth = if (is.null(params$linewidth)) 1 else params$linewidth,
  )
}



geom_are_subclones <- function(resid, params) {
  dt <- resid |>
    pivot_longer(
      cols = c("Clone", starts_with("Subclone")),
      names_to = "component",
      values_to = "pred"
    ) |>
    filter(!is.na(.data$pred))
  geom_area(
    aes(.data$f, .data$pred, group = .data$component),
    data = dt,
    position = "identity",
    alpha = if (is.null(params$alpha)) 0.3 else params$alpha,
    color = if (is.null(params$color)) "black" else params$color,
    linetype = if (is.null(params$linetype)) "solid" else params$linetype,
    linewidth = if (is.null(params$linewidth)) 1 else params$linewidth,
    show.legend = FALSE
  )
}



geom_line_final_fit <- function(resid, params) {
  geom_line(
    aes(.data$f, .data$model_pred),
    data = resid |> filter(.data$powerlaw_pred < .data$ylim),
    alpha = if (is.null(params$alpha)) 1 else params$alpha,
    color = if (is.null(params$color)) "red" else params$color,
    linetype = if (is.null(params$linetype)) "solid" else params$linetype,
    linewidth = if (is.null(params$linewidth)) 1 else params$linewidth
  )
}


## --------------------------- Other functions --------------------------------

plot_clones <- function(clones) {
  clones |>
    get_binomial_predictions() |>
    select(-"binom_pred") |>
    pivot_longer(-.data$f) |>
    ggplot(aes(.data$f, .data$value, group = .data$name)) +
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
    join_aes(aes(.data$f, .data$y), mapping),
    data = . %>%
      mutate(y = A * 1 / .data$f ^ alpha) %>%
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
    map(~get_model_residuals(object, .x)) |>
    bind_rows(.id = "model_name") |>
    left_join(ylimits, by = "sample_id") |>
    join_metadata(object) |>
    filter(!!sym(column_name) < .data$ylim, .data$f >= 0)

  plot_SFS(object, geom = "bar", ...) +
    geom_line(
      aes(y = !!sym(column_name), color = .data$model_name, group = .data$model_name),
      data = resids,
      linetype = linetype, linewidth = linewidth
    ) +
    facet_wrap(~.data$sample_id, scales = "free_y")
}

