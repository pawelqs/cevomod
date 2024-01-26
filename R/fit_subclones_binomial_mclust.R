

#' @describeIn fit_subclones Fit subclonal distributions to neutral model residuals using mclust
#' @export
fit_subclones_mclust <- function(object,
                                 N = 1:3,
                                 powerlaw_model_name = active_models(object),
                                 snvs_name = default_SNVs(object),
                                 upper_f_limit = 0.75,
                                 verbose = get_verbosity()) {
  msg("Fitting binomial models using mclust", verbose = verbose)

  powerlaw_models <- get_models(object, powerlaw_model_name)
  stop_if_models_not_powerlaw(powerlaw_models, powerlaw_model_name)

  residuals <- get_model_residuals(object, model_name = powerlaw_model_name) |>
    filter(.data$f >= 0)
  sequencing_depths <- SNVs(object, name = snvs_name) |>
    get_local_sequencing_depths() |>
    transmute(.data$sample_id, .data$f, sequencing_DP = .data$median_DP)

  pb <- if (verbose) progress_bar$new(total = n_distinct(residuals$sample_id)) else NULL
  coefs <- residuals |>
    select("sample_id", "f", "powerlaw_resid_clones") |>
    mutate(
      powerlaw_resid_clones = if_else(.data$f > upper_f_limit, 0, .data$powerlaw_resid_clones)
    ) |>
    nest_by(.data$sample_id) |>
    reframe(fit_binomial_models(.data$data, N = N, pb = pb)) |>
    mutate(model = "binomial_clones", .after = "sample_id") |>
    mutate(f = round(.data$frequency, digits = 2)) |>
    left_join(sequencing_depths, by = c("sample_id", "f")) |>
    select(-"f") |>
    evaluate_binomial_models()

  coefs
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
  clones <- tibble(
    N = length(mclust_model$parameters$mean),
    frequency = mclust_model$parameters$mean,
    N_mutations = round(mclust_model$parameters$pro * n_mutations),
    BIC = mclust_model$bic
  ) |>
    arrange(desc(.data$frequency)) |>
    mutate(
      component = if_else(row_number() == 1, "Clone", str_c("Subclone ", row_number() - 1)),
      .before = "frequency"
    )
  names(clones$frequency) <- NULL
  names(clones$BIC) <- NULL
  clones
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

