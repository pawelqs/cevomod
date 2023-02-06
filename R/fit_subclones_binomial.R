
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
#' @param powerlaw_model residual of which powerlaw model to use?
#'   williams_neutral/tung_durrett
#' @param verbose verbose?
#' @export
fit_subclones.cevodata <- function(object,
                                   N = 1:3,
                                   powerlaw_model_name = active_models(object),
                                   verbose = TRUE, ...) {
  msg("Fitting binomial models", verbose = verbose)
  powerlaw_models <- get_powerlaw_models(object, powerlaw_model_name)

  residuals <- get_residuals(object, model = powerlaw_model_name) |>
    filter(.data$VAF >= 0 )

  sequencing_depths <- SNVs(object) |>
    get_local_sequencing_depths() |>
    transmute(.data$sample_id, .data$VAF, sequencing_DP = .data$median_DP)

  models <- residuals |>
    select("sample_id", "VAF", "powerlaw_resid_clones") |>
    nest_by(.data$sample_id) |>
    summarise(fit_binomial_models(.data$data, N = N), .groups = "drop") |>
    mutate(model = "binomial_clones", .after = "sample_id") |>
    mutate(VAF = round(.data$cellularity, digits = 2)) |>
    left_join(sequencing_depths, by = c("sample_id", "VAF")) |>
    select(-"VAF") |>
    evaluate_binomial_models()

  best_models <- models |>
    filter(.data$best) |>
    nest_by(.data$sample_id, .key = "clones")

  clonal_predictions <- residuals |>
    select("sample_id", "VAF_interval", "VAF") |>
    nest_by(.data$sample_id, .key = "VAFs") |>
    left_join(best_models, by = "sample_id") |>
    summarise(get_binomial_predictions(.data$clones, .data$VAFs), .groups = "drop") |>
    select(-"VAF")

  residuals <- residuals |>
    select(-"model_resid") |>
    left_join(clonal_predictions, by = c("sample_id", "VAF_interval")) |>
    mutate(
      model_pred = .data$powerlaw_pred + .data$binom_pred,
      model_resid = .data$SFS - .data$model_pred
    )

  models <- powerlaw_models |>
    bind_rows(models) |>
    arrange(.data$sample_id, .data$best, .data$model)

  models_name <- paste0(powerlaw_model_name, "_subclones")
  resid_name <- paste0("residuals_", models_name)
  object$models[[models_name]] <- models
  object$misc[[resid_name]] <- residuals
  object$active_models <- models_name
  object
}


fit_binomial_models <- function(...) {
  fit_binomial_models_Mclust(...)
}


fit_binomial_models_Mclust <- function(residuals, N) {
  VAFs <- rep(residuals$VAF, times = floor(residuals$powerlaw_resid_clones))
  if (length(VAFs) == 0) {
    return(empty_clones_tibble())
  }

  mclust_res <- N |>
    map(~ mclust::Mclust(VAFs, G = .x, verbose = FALSE)) |>
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
    unnest("data") |>
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
    pmap(get_binomial_distribution) |>
    map(rebinarize_distribution, n_bins = 100) |>
    map("pred")
  names(x) <- clones$component
  x <- as_tibble(x) |>
    corrr::correlate(quiet = TRUE) |>
    select(-"term")
  x[is.na(x)] <- 0
  any(x > 0.5)
}


get_binomial_predictions <- function(clones, VAFs) {
  clones_predictions <- clones |>
    pmap(get_binomial_distribution) |>
    map(rebinarize_distribution, VAFs = VAFs$VAF) |>
    map("pred") |>
    set_names(clones$component) |>
    bind_cols()
  res <- bind_cols(
    VAFs,
    clones_predictions,
    binom_pred = rowSums(clones_predictions)
  )
}


get_binomial_distribution <- function(cellularity, N_mutations, sequencing_DP, ...) {
  i <- 0:round(sequencing_DP)
  tibble(
    VAF = i / sequencing_DP,
    pred = N_mutations * stats::dbinom(i, round(sequencing_DP), cellularity)
  )
}


rebinarize_distribution <- function(distribution, n_bins = NULL, VAFs = NULL) {
  if (is.null(n_bins) == is.null(VAFs)) {
    stop("Provide n_bins OR VAFs")
  }
  new_VAFs <- if (is.null(VAFs)) (1:n_bins) / n_bins else VAFs
  original_VAFs <- distribution$VAF
  distribution$VAF <- NULL

  new_distributions <- distribution |>
    map_dfc(~ stats::approx(original_VAFs, .x, xout = new_VAFs, rule = 2)$y)
  scaling_factors <- map2(new_distributions, distribution, ~ sum(.x) / sum(.y))
  rescaled_distributions <- map2_dfc(new_distributions, scaling_factors, ~ .x / .y)

  rescaled_distributions |>
    mutate(VAF = new_VAFs) |>
    select("VAF", everything())
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
