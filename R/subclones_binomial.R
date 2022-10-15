
#' Fit subclonal distributions to neutral model residuals
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
  # TODO: smooth(meutral_resid)?

  clones <- residuals |>
    group_by(.data$sample_id) |>
    nest() |>
    mutate(clones = map(.data$data, fit_binomial_models_Mclust, clones = 1:3, epochs = 100, eps = 1e-3)) |>
    select(-.data$data) |>
    unnest(.data$clones) |>
    ungroup()

  clonal_predictions <- clones |>
    nest_by(sample_id) |>
    deframe() |>
    map(get_binomial_predictions) |>
    bind_rows(.id = "sample_id")

  residuals <- object$models$residuals |>
    left_join(clonal_predictions, by = c("sample_id", "VAF")) |>
    mutate(
      model_pred = .data$neutral_pred + .data$binom_pred,
      model_resid = .data$SFS - .data$model_pred
    )

  object$models[["binomial_models"]] <- clones
  object$models[["residuals"]] <- residuals
  object$SNVs[[default_SNVs(object)]] <- classify_SNVs(SNVs(object), residuals)
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
      BIC = model$bic,
      best = TRUE
    ) |>
    arrange(desc(cellularity)) |>
    mutate(
      clone = if_else(row_number() == 1, "Clone", str_c("Subclone ", row_number() - 1)),
      .before = "cellularity"
    )
  clones
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


get_binomial_predictions <- function(clones) {
  predictions <- tibble(
    i = 1:100,
    VAF = .data$i/100,
    subclonal_pred = map2(clones$N_mutations, clones$cellularity, ~.x * dbinom(.data$i, 100, .y)) |>
      set_names(clones$clone) |>
      as_tibble(),
    binom_pred = rowSums(subclonal_pred)
  )
  predictions |>
    unnest(subclonal_pred) |>
    select(-.data$i)
}


# predict_binomial_distribution <- function(Ns, means) {
#   tibble(
#     i = 1:100,
#     VAF = .data$i/100,
#     binom_pred = map2(unlist(Ns), unlist(means), ~.x * dbinom(.data$i, 100, .y)) |>
#       reduce(`+`)
#   )
# }


classify_SNVs <- function(snvs, residuals) {
  probabilities <- get_probabilities_tbl(residuals)
  snvs |>
    mutate(VAF_chr = as.character(round(VAF, digits = 2))) |>
    left_join(probabilities, by = c("sample_id", "VAF_chr")) |>
    select(-.data$VAF_chr)
}


get_probabilities_tbl <- function(residuals) {
  probabilities <- residuals |>
    mutate(VAF = as.character(VAF)) |>
    select(sample_id, VAF, Neutral = neutral_pred, Clone, starts_with("Subclone"), model_pred) |>
    transmute(
      sample_id,
      VAF_chr = VAF,
      across(c("Neutral", "Clone", starts_with("Subclone")), ~.x/model_pred)
    )
  probabilities
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
                                 binomial_layer = FALSE,
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
    mutate(ylim = max(.data$SFS) * 1.2) |>
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
    if (binomial_layer && subclones_fitted) {
      geom_line(
        aes(.data$VAF, .data$binom_pred),
        data = resid,
        size = 1, color = "black"
      )
    },
    if (subclones && subclones_fitted) {
      dt <- resid |>
        pivot_longer(
          cols = c("Clone", starts_with("Subclone")),
          names_to = "clone",
          values_to = "pred"
        ) |>
        filter(!is.na(pred))
      geom_line(
        aes(.data$VAF, .data$pred, group = clone),
        data = dt,
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
