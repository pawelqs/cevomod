
#' Fit subclonal distributions to neutral model residuals
#'
#' @param object object
#' @param ... other arguments
#' @export
fit_subclones <- function(object, ...) {
  UseMethod("fit_subclones")
}


#' @export
fit_subclones.cevodata <- function(object, N = 1:3, ...) {
  residuals <- get_residuals(object, model = "neutral_models")

  clones <- residuals |>
    nest_by(.data$sample_id) |>
    summarise(
      model = "binomial_clones",
      clones = fit_binomial_models_Mclust(.data$data, N = N),
      .groups = "drop"
    ) |>
    unnest(.data$clones) |>
    mutate(VAF = round(.data$cellularity, digits = 2)) |>
    left_join(get_sequencing_depths(object), by = c("sample_id", "VAF")) |>
    select(-.data$VAF)

  clonal_predictions <- clones |>
    nest_by(.data$sample_id) |>
    deframe() |>
    map(get_binomial_predictions) |>
    bind_rows(.id = "sample_id")

  residuals <- residuals |>
    left_join(clonal_predictions, by = c("sample_id", "VAF")) |>
    mutate(
      model_pred = .data$neutral_pred + .data$binom_pred,
      model_resid = .data$SFS - .data$model_pred
    )

  models <- bind_rows(
    get_neutral_models(object),
    clones
  ) |>
    arrange(.data$sample_id)

  object$models[["binomial_models"]] <- clones
  object$residuals[["binomial_models"]] <- residuals
  # object$SNVs[[default_SNVs(object)]] <- classify_SNVs(SNVs(object), residuals)
  object$active_model <- "binomial_models"
  object
}


fit_binomial_models_Mclust <- function(residuals, N) {
  VAFs <- rep(residuals$VAF, times = floor(residuals$neutral_resid_clones))
  clones <- N |>
    map(~mclust::Mclust(VAFs, G = .x, verbose = FALSE)) |>
    discard(any_clusters_overlap) |>
    map(mclust_to_clones_tbl, n_mutations = length(VAFs)) |>
    bind_rows() |>
    mutate(best = .data$BIC == max(.data$BIC))
  clones |>
    filter(.data$best)
}


any_clusters_overlap <- function(mclust_res) {
  mean <- mclust_res$parameters$mean
  sd <- sqrt(mclust_res$parameters$variance$sigmasq)
  clust_ranges <- IRanges::IRanges(
    start = round(100 * (mean - sd)),
    end = round(100 * (mean + sd)),
  )
  range_coverage <- IRanges::coverage(clust_ranges)
  any(range_coverage@values > 1)
}


mclust_to_clones_tbl <- function(mclust_model, n_mutations) {
  tibble(
    n = length(mclust_model$parameters$mean),
    cellularity = mclust_model$parameters$mean,
    N_mutations = round(mclust_model$parameters$pro * n_mutations),
    BIC = mclust_model$bic
  ) |>
    arrange(desc(.data$cellularity)) |>
    mutate(
      component = if_else(row_number() == 1, "Clone", str_c("Subclone ", row_number() - 1)),
      .before = "cellularity"
    )
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
  clones_predictions <- clones |>
    rename(sequencing_DP = .data$median_DP) |>
    pmap(get_binomial_distribution) |>
    map(rebinarize_distribution, n_bins = 100) |>
    set_names(clones$component) |>
    bind_rows(.id = "component") |>
    pivot_wider(names_from = "component", values_from = "pred")
  clones_predictions$binom_pred <- clones_predictions |>
    select(-.data$VAF) |>
    rowSums()
  clones_predictions
  # predictions <- tibble(
  #   i = 1:100,
  #   VAF = .data$i/100,
  #   subclonal_pred = map2(clones$N_mutations, clones$cellularity, ~.x * dbinom(.data$i, 100, .y)) |>
  #     set_names(clones$component) |>
  #     as_tibble(),
  #   binom_pred = rowSums(.data$subclonal_pred)
  # )
  # predictions |>
  #   unnest(.data$subclonal_pred) |>
  #   select(-.data$i)
}


get_binomial_distribution <- function(cellularity, N_mutations, sequencing_DP, ...) {
  i <- 0:round(sequencing_DP)
  tibble(
    VAF = i/sequencing_DP,
    pred = N_mutations * stats::dbinom(i, round(sequencing_DP), cellularity)
  )
}


rebinarize_distribution <- function(distribution, n_bins = 100) {
  i <- 1:n_bins
  new_distribution <- tibble(
    VAF = i/n_bins,
    pred = stats::approx(distribution$VAF, distribution$pred, xout = .data$VAF, rule = 2)$y
  )
  scaling_factor <- sum(new_distribution$pred) / sum(distribution$pred)
  new_distribution |>
    mutate(pred = .data$pred / scaling_factor)
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
    mutate(VAF_chr = as.character(round(.data$VAF, digits = 2))) |>
    left_join(probabilities, by = c("sample_id", "VAF_chr")) |>
    select(-.data$VAF_chr)
}


get_probabilities_tbl <- function(residuals) {
  probabilities <- residuals |>
    mutate(VAF = as.character(.data$VAF)) |>
    select(
      .data$sample_id, .data$VAF,
      Neutral = .data$neutral_pred,
      .data$Clone,
      starts_with("Subclone"),
      .data$model_pred
    ) |>
    transmute(
      .data$sample_id,
      VAF_chr = .data$VAF,
      across(c("Neutral", "Clone", starts_with("Subclone")), ~.x/model_pred)
    )
  probabilities
}


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
    select(.data$sample_id, .data$from, .data$to)

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
        data = resid, # |> filter(.data$neutral_pred < .data$ylim),
         fill = "white", color = "white",
        alpha = 0.3,
        size = 0.5, show.legend = FALSE
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
        size = 1, color = "black"
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
