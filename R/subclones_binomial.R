
#' fit_subclones
#'
#' @param object object
#' @param ... other arguments
#' @export
fit_subclones <- function(object, ...) {
  UseMethod("fit_subclones")
}


#' @export
fit_subclones.cevodata <- function(object, ...) {
  residuals <- object$models$residuals
  residuals$binom_pred <- NULL
  #   mutate(sm = smooth(meutral_resid))

  models <- residuals |>
    group_by(.data$sample_id) |>
    nest() |>
    mutate(clones = map(.data$data, fit_binomial_models_Mclust, clones = 1:3, epochs = 100, eps = 1e-3)) |>
    select(-.data$data) |>
    unnest(.data$clones)
  # best_model <- models |>
  #   slice(1)

  binom_residuals <- models |>
    group_by(sample_id) |>
    summarise(prediction = predict_binoms(N_mutations, cellularity), .groups = "drop") |>
    unnest(.data$prediction) |>
    mutate(VAF = as.character(.data$VAF)) |>
    select(-.data$i)
  residuals <- residuals |>
    mutate(VAF = as.character(.data$VAF)) |>
    left_join(binom_residuals, by = c("sample_id", "VAF")) |>
    mutate(
      VAF = parse_double(.data$VAF),
      model_pred = .data$neutral_pred + .data$binom_pred,
      model_resid = .data$y - .data$model_pred
    )

  object$models[["binomial_models"]] <- models
  object$models[["residuals"]] <- residuals
  object
}


fit_binomial_models_Mclust <- function(residuals, clones, epochs, eps) {
  data <- residuals |>
    mutate() |>
    transmute(
      .data$VAF,
      n = as.integer(round(.data$neutral_resid_clones)),
      muts = map(.data$n, ~tibble(i = 1:.x))
    ) |>
    unnest(.data$muts)
  model <- mclust::Mclust(data$VAF, G = 1:3, verbose = FALSE)
  clones <- tibble(
      cellularity = model$parameters$mean,
      N_mutations = round(model$parameters$pro * nrow(data)),
      clone = if_else(cellularity == max(cellularity), "Clone", str_c("Subclone ", row_number())),
      BIC = model$bic,
      best = TRUE
    )
  clones |>
    select(clone, everything())
}


# fit_binomial_models <- function(er, clones, epochs, eps) {
#   clones |>
#     set_names(clones) |>
#     map(~.fit_binomial_models(er, .x, epochs, eps)) |>
#     bind_rows(.id = "clones")
# }
#
#
# .fit_binomial_models <- function(er, clones, epochs, eps) {
#   tibble(
#     means = list(runif(clones)),
#     Ns = 300,
#     BIC = runif(1)
#   )
# }


predict_binoms <- function(Ns, means) {
  tibble(
    i = 1:100,
    VAF = .data$i/100,
    binom_pred = map2(unlist(Ns), unlist(means), ~.x * dbinom(.data$i, 100, .y)) |>
      reduce(`+`)
  )
}


#' Plot cevodata models
#' @param object cevodata object
#' @param neutral_tail TRUE,
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
                                 subclones = TRUE,
                                 final_fit = TRUE,
                                 ...) {

  neutral_lm_fitted <- !is.null(object$models$neutral_lm)
  subclones_fitted <- !is.null(object$models$residuals$binom_pred)

  lm_models <- object$models$neutral_lm |>
    filter(.data$best) |>
    select(.data$sample_id, .data$from, .data$to)

  resid <- object$models$residuals |>
    left_join(lm_models, by = "sample_id") |>
    group_by(.data$sample_id) |>
    mutate(ylim = max(.data$y) * 1.2) |>
    ungroup() |>
    mutate(neutr = (.data$VAF >= .data$from & .data$VAF <= .data$to))

  model_layers <- list(
    if (neutral_tail && neutral_lm_fitted) {
      geom_line(
        aes(.data$VAF, .data$neutral_pred),
        data = resid |> filter(.data$neutral_pred < .data$ylim),
        color = "black", size = 1, linetype = "dashed", show.legend = FALSE
      )
    },
    if (neutral_tail && neutral_lm_fitted) {
      geom_line(
        aes(.data$VAF, .data$neutral_pred),
        data = resid |> filter(.data$neutr, .data$neutral_pred < .data$ylim),
        color = "black", size = 1, show.legend = FALSE
      )
    },
    if (subclones && subclones_fitted) {
      geom_line(
        aes(.data$VAF, .data$binom_pred),
        data = resid,
        size = 1, color = "black"
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
