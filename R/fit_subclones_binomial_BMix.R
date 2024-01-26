
#' @describeIn fit_subclones Fit subclonal distributions to neutral model residuals using BMix
#' @export
fit_subclones_bmix <- function(object,
                               N = 1:3,
                               powerlaw_model_name = active_models(object),
                               snvs_name = default_SNVs(object),
                               upper_f_limit = 0.75,
                               verbose = get_verbosity()) {
  msg("Fitting binomial models using BMix", verbose = verbose)
  rlang::check_installed("BMix", reason = "to fit subclomes using BMix")

  powerlaw_models <- get_models(object, powerlaw_model_name)
  stop_if_models_not_powerlaw(powerlaw_models, powerlaw_model_name)

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

  coefs
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

