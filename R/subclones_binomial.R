
#' Fit subclonal distributions to neutral model residuals
#'
#' @param object object
#' @param ... other arguments
#' @export
fit_subclones <- function(object, ...) {
  UseMethod("fit_subclones")
}


#' @rdname fit_subclones
#' @param N vector of numbers of clones to test
#' @export
fit_subclones.cevodata <- function(object, N = 1:3, ...) {
  residuals <- get_residuals(object, model = "neutral_models")

  models <- residuals |>
    nest_by(.data$sample_id) |>
    summarise(fit_binomial_models(.data$data, N = N), .groups = "drop") |>
    mutate(model = "binomial_clones", .after = "sample_id") |>
    mutate(VAF = round(.data$cellularity, digits = 2)) |>
    left_join(get_sequencing_depths(object), by = c("sample_id", "VAF")) |>
    select(-.data$VAF) |>
    evaluate_binomial_models()

  clonal_predictions <- models |>
    filter(.data$best) |>
    nest_by(.data$sample_id) |>
    summarise(get_binomial_predictions(.data$data), .groups = "drop")

  residuals <- residuals |>
    left_join(clonal_predictions, by = c("sample_id", "VAF")) |>
    mutate(
      model_pred = .data$neutral_pred + .data$binom_pred,
      model_resid = .data$SFS - .data$model_pred
    )

  models <- get_neutral_models(object) |>
    bind_rows(models) |>
    arrange(.data$sample_id, .data$best, .data$model)

  object$models[["binomial_models"]] <- models
  object$residuals[["binomial_models"]] <- residuals
  # object$SNVs[[default_SNVs(object)]] <- classify_SNVs(SNVs(object), residuals)
  object$active_models <- "binomial_models"
  object
}


fit_binomial_models <- function(...) {
  fit_binomial_models_Mclust(...)
}


fit_binomial_models_Mclust <- function(residuals, N) {
  VAFs <- rep(residuals$VAF, times = floor(residuals$neutral_resid_clones))
  if (length(VAFs) == 0) {
    return(empty_clones_tibble())
  }

  mclust_res <- N |>
    map(~mclust::Mclust(VAFs, G = .x, verbose = FALSE)) |>
    discard(is.null)

  if (length(mclust_res) > 0) {
    clones <- mclust_res |>
      map(mclust_to_clones_tbl, n_mutations = length(VAFs)) |>
      bind_rows()
  } else {
    clones <- empty_clones_tibble()
  }
  clones
}


mclust_to_clones_tbl <- function(mclust_model, n_mutations) {
  tibble(
    N = length(mclust_model$parameters$mean),
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


empty_clones_tibble <- function() {
  tibble(
    N = integer(),
    component = character(),
    cellularity = double(),
    N_mutations = double(),
    BIC = double()
  )
}


evaluate_binomial_models <- function(models) {
  models <- models |>
    nest_by(.data$sample_id, .data$model, .data$N) |>
    mutate(has_overlapping_clones = any_binomial_distibutions_correlate(.data$data)) |>
    unnest(.data$data) |>
    ungroup()
  best_BIC_values <- models |>
    filter(!.data$has_overlapping_clones) |>
    group_by(.data$sample_id) |>
    summarise(max_BIC = max(.data$BIC))
  models |>
    left_join(best_BIC_values, by = "sample_id") |>
    mutate(best = (.data$BIC == .data$max_BIC) & !.data$has_overlapping_clones)
}


any_binomial_distibutions_correlate <- function(clones) {
  if (nrow(clones) == 1) {
    return(FALSE)
  }
  x <- clones |>
    rename(sequencing_DP = .data$median_DP) |>
    pmap(get_binomial_distribution) |>
    map(rebinarize_distribution, n_bins = 100) |>
    map("pred")
  names(x) <- clones$component
  x <- as_tibble(x) |>
    corrr::correlate(quiet = TRUE) |>
    select(-.data$term)
  x[is.na(x)] <- 0
  any(x > 0.5)
}


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


were_subclonal_models_fitted <- function(object, ...) {
  models <- get_models(object)
  expect_colnames <- c("N", "cellularity", "N_mutations")
  all(expect_colnames %in% colnames(models))
}
