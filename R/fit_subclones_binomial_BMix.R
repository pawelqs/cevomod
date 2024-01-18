
#' @describeIn fit_subclones Fit subclonal distributions to neutral model residuals using BMix
#' @export
fit_subclones_bmix <- function(object,
                               N = 1:3,
                               powerlaw_model_name = active_models(object),
                               name = paste0(powerlaw_model_name, "_subclones"),
                               snvs_name = default_SNVs(object),
                               upper_f_limit = 0.75,
                               verbose = get_verbosity()) {
  msg("Fitting binomial models using BMix", verbose = verbose)
  rlang::check_installed("BMix", reason = "to fit subclomes using BMix")

  powerlaw_models <- get_models(object, powerlaw_model_name)
  if ("cv_powerlaw_models" %not in% class(powerlaw_models)) {
    stop(
      powerlaw_model_name, " is not a powerlaw model, required to fit subclones.",
      "Use another model"
    )
  }
  residuals <- get_model_residuals(object, model_name = powerlaw_model_name) |>
    filter(.data$f >= 0)
  # necessary, since BMix does not return this parameter
  sequencing_depth <- SNVs(object, name = snvs_name) |>
    get_local_sequencing_depths() |>
    transmute(.data$sample_id, .data$f, sequencing_DP = .data$median_DP)

  non_neutral_tail_mut_counts <- residuals |>
    mutate(n = if_else(.data$f > upper_f_limit, 0, round(.data$powerlaw_resid_clones))) |>
    select("sample_id", "f_interval", "n")
  snvs_to_cluster <- SNVs(object, name = snvs_name) |>
    nest_by(.data$sample_id, .data$f_interval) |>
    inner_join(non_neutral_tail_mut_counts, by = c("sample_id", "f_interval")) |>
    reframe(
      slice_sample(.data$data, n = .data$n)
    )

  data <- snvs_to_cluster |>
    select("sample_id", "alt_reads", "DP") |>
    nest_by(.data$sample_id) |>
    mutate(data = list(as.data.frame(.data$data)))
  pb <- if (verbose) progress_bar$new(total = nrow(data)) else NULL
  coefs <- data |>
    reframe(fit_binomial_models_BMix(.data$data, N, pb, verbose)) |>
    mutate(model = "binomial_clones_BMix", .after = "sample_id") |>
    mutate(
      f = round(.data$cellularity, digits = 2),
      best = TRUE
    ) |>
    left_join(sequencing_depth, by = c("sample_id", "f")) |>
    select(-"f")

  best_model_coefs <- coefs |>
    filter(.data$best) |>
    nest_by(.data$sample_id, .key = "clones")

  clonal_predictions <- residuals |>
    select("sample_id", "f_interval", "f") |>
    nest_by(.data$sample_id, .key = "intervals") |>
    inner_join(best_model_coefs, by = "sample_id") |>
    reframe(get_binomial_predictions(.data$clones, .data$intervals)) |>
    select(-"f")

  residuals <- residuals |>
    select(-"model_resid") |>
    left_join(clonal_predictions, by = c("sample_id", "f_interval")) |>
    mutate(
      model_pred = .data$powerlaw_pred + .data$binom_pred,
      model_resid = .data$SFS - .data$model_pred
    )

  coefs <- powerlaw_models$coefs |>
    bind_rows(coefs) |>
    arrange(.data$sample_id, .data$best, .data$model)

  models <- lst(coefs, residuals, info = powerlaw_models$info)
  # models_name <- paste0(powerlaw_model_name, "_subclones")
  # resid_name <- paste0("residuals_", models_name)
  # object$models[[models_name]] <- models
  # object$misc[[resid_name]] <- residuals
  # object$active_models <- models_name
  # object
  add_models(object, models, name)
}


fit_binomial_models_BMix <- function(data, N, pb = NULL, verbose = TRUE) {
  bmixfit <- safely(BMix::bmixfit)
  model <- bmixfit(data, K.Binomials = N, K.BetaBinomials = 0, silent = verbose < 2)

  if (is.null(model$result)) {
    model_tbl <- empty_clones_tibble()
  } else {
    model <- model$result
    mut_counts <- tibble(cluster = model$labels) |>
      count(.data$cluster, name = "N_mutations")
    model_tbl <- BMix::Parameters(model) |>
      arrange(desc(.data$mean)) |>
      mutate(
        N = n(),
        component = if_else(row_number() == 1, "Clone", str_c("Subclone ", row_number() - 1)),
        BIC = model$BIC
      ) |>
      left_join(mut_counts, by = "cluster") |>
      select("N", "component", cellularity = "mean", "N_mutations", "BIC")
  }

  if (!is.null(pb)) pb$tick()
  model_tbl
}

