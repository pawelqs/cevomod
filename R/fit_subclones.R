

#' Fit clonal and subclonal components of the model to the residuals of the
#' power-law model
#'
#' @param object cevodata object
#' @param N Vector of numbers of clones to for models
#' @param powerlaw_model_name Residual of which powerlaw model to use?
#'   powerlaw_fixed/powerlaw_optim
#' @param snvs_name Which snvs to to use?
#' @param cnas_name Which cnas to to use?
#' @param method Clustering method to use. Currently supported methods:
#'   - mclust - the fastest method, approximately 3-4 times faster than BMix,
#'     but uses a gaussian mixture modelling
#'   - BMix - is more accurate, considers subclones as binomial clusters,
#'     slightly slower
#'   - CliP - Clonal structure identification through penalizing pairwise
#'     differences
#' @param upper_f_limit ignore variants with f higher than
#' @param verbose Verbose?
#'
#' @examples
#' \dontrun{
#' # Using BMix
#' fit_subclones(test_data_fitted)
#' # or
#' fit_subclones_bmix(test_data_fitted)
#'
#' # Using mclust
#' fit_subclones(test_data_fitted, method = "mclust")
#' # or
#' fit_subclones_mclust(test_data_fitted)
#'
#' # Using CliP
#' set_containers_dir(selected_dir)
#' build_clip_container()
#' fit_subclones(test_data_fitted, method = "CliP")
#' # or
#' fit_subclones_clip(test_data_fitted)
#' }
#' @name fit_subclones
NULL



#' @describeIn fit_subclones Provides a common interface for all other methods,
#'   runs the selected method and passes all the required arguments down.
#' @export
fit_subclones <- function(object,
                          N = 1:3,
                          powerlaw_model_name = active_models(object),
                          name = paste0(powerlaw_model_name, "_subclones"),
                          snvs_name = default_SNVs(object),
                          cnas_name = default_CNAs(object),
                          method = "BMix",
                          upper_f_limit = 0.75,
                          clip_sif = NULL,
                          clip_input = file.path(tempdir(), "clip_input"),
                          clip_output = file.path(tempdir(), "clip_output"),
                          verbose = get_verbosity()) {
  powerlaw_models <- get_models(object, powerlaw_model_name)
  stop_if_models_not_powerlaw(powerlaw_models, powerlaw_model_name)

  residuals <- get_model_residuals(object, model_name = powerlaw_model_name) |>
    filter(.data$f >= 0)

  if (method == "BMix") {
    coefs <- object |>
      fit_subclones_bmix(N, powerlaw_model_name, snvs_name, upper_f_limit, verbose)
  } else if (method == "mclust") {
    coefs <- object |>
      fit_subclones_mclust(N, powerlaw_model_name, snvs_name, upper_f_limit, verbose)
  } else if (method == "CliP") {
    coefs <- object |>
      fit_subclones_clip(powerlaw_model_name, snvs_name, cnas_name, upper_f_limit, verbose = verbose)
  } else {
    stop("Currently supported methods are: BMix, CliP, and mclust")
  }

  best_coefs <- coefs |>
    filter(.data$best) |>
    nest_by(.data$sample_id, .key = "clones")

  clonal_predictions <- residuals |>
    select("sample_id", "f_interval", "f") |>
    nest_by(.data$sample_id, .key = "intervals") |>
    inner_join(best_coefs, by = "sample_id") |>
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
  class(models) <- c("cv_powerlaw_subclones_models", "list")
  add_models(object, models, name)
}


get_binomial_predictions <- function(clones, intervals) {
  clones_predictions <- clones |>
    pmap(get_binomial_distribution) |>
    map(rebinarize_distribution, f = intervals$f) |>
    map("pred") |>
    set_names(clones$component) |>
    bind_cols()
  res <- bind_cols(
    intervals,
    clones_predictions,
    binom_pred = rowSums(clones_predictions)
  )
}


get_binomial_distribution <- function(cellularity, N_mutations, sequencing_DP, ...) {
  i <- 0:round(sequencing_DP)
  tibble(
    f = i / sequencing_DP,
    pred = N_mutations * stats::dbinom(i, round(sequencing_DP), cellularity)
  )
}


rebinarize_distribution <- function(distribution, n_bins = NULL, f = NULL) {
  if (is.null(n_bins) == is.null(f)) {
    stop("Provide n_bins OR f")
  }
  new_fs <- if (is.null(f)) (1:n_bins) / n_bins else f
  original_fs <- distribution$f
  distribution$f <- NULL

  new_distributions <- distribution |>
    map_dfc(~ stats::approx(original_fs, .x, xout = new_fs, rule = 2)$y)
  scaling_factors <- map2(new_distributions, distribution, ~ sum(.x) / sum(.y))
  rescaled_distributions <- map2_dfc(new_distributions, scaling_factors, ~ .x / .y)

  rescaled_distributions |>
    mutate(f = new_fs) |>
    select("f", everything())
}


were_subclonal_models_fitted <- function(object, ...) {
  models <- get_model_coefficients(object)
  expect_colnames <- c("N", "cellularity", "N_mutations")
  all(expect_colnames %in% colnames(models))
}


stop_if_models_not_powerlaw <- function(models, name) {
  if ("cv_powerlaw_models" %not in% class(models)) {
    stop(
      name, " is not a powerlaw model, required to fit subclones.",
      "Use another model"
    )
  }
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
