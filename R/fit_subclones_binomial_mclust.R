

#' @describeIn fit_subclones Fit subclonal distributions to neutral model residuals using mclust
#' @export
fit_subclones_mclust <- function(object,
                                 N = 1:3,
                                 powerlaw_model_name = active_models(object),
                                 snvs_name = default_SNVs(object),
                                 upper_f_limit = 0.75,
                                 verbose = get_cevomod_verbosity()) {
  msg("Fitting binomial models using mclust", verbose = verbose)

  powerlaw_models <- get_powerlaw_models(object, powerlaw_model_name)
  residuals <- get_residuals(object, models_name = powerlaw_model_name) |>
    filter(.data$f >= 0)
  sequencing_depths <- SNVs(object, which = snvs_name) |>
    get_local_sequencing_depths() |>
    transmute(.data$sample_id, .data$f, sequencing_DP = .data$median_DP)

  pb <- if (verbose) progress_bar$new(total = n_distinct(residuals$sample_id)) else NULL
  models <- residuals |>
    select("sample_id", "f", "powerlaw_resid_clones") |>
    mutate(
      powerlaw_resid_clones = if_else(.data$f > upper_f_limit, 0, .data$powerlaw_resid_clones)
    ) |>
    nest_by(.data$sample_id) |>
    reframe(fit_binomial_models(.data$data, N = N, pb = pb)) |>
    mutate(model = "binomial_clones", .after = "sample_id") |>
    mutate(f = round(.data$cellularity, digits = 2)) |>
    left_join(sequencing_depths, by = c("sample_id", "f")) |>
    select(-"f") |>
    evaluate_binomial_models()

  best_models <- models |>
    filter(.data$best) |>
    nest_by(.data$sample_id, .key = "clones")

  clonal_predictions <- residuals |>
    select("sample_id", "f_interval", "f") |>
    nest_by(.data$sample_id, .key = "intervals") |>
    inner_join(best_models, by = "sample_id") |>
    reframe(get_binomial_predictions(.data$clones, .data$intervals)) |>
    select(-"f")

  residuals <- residuals |>
    select(-"model_resid") |>
    left_join(clonal_predictions, by = c("sample_id", "f_interval")) |>
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


fit_binomial_models <- function(..., pb = NULL) {
  fit <- fit_binomial_models_Mclust(...)
  if (!is.null(pb)) pb$tick()
  fit
}


fit_binomial_models_Mclust <- function(residuals, N) {
  f <- rep(residuals$f, times = floor(residuals$powerlaw_resid_clones))
  if (length(f) <= 1) {
    return(empty_clones_tibble())
  }
  if (sd(f) == 0) {
    N <- 1
  }
  mclust_res <- N |>
    map(~ mclust::Mclust(f, G = .x, verbose = FALSE)) |>
    discard(is.null)

  if (length(mclust_res) > 0) {
    clones <- mclust_res |>
      map(mclust_to_clones_tbl, n_mutations = length(f)) |>
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

