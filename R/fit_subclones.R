

#' Fit clonal and subclonal components of the model to the residuals of the
#' power-law model
#'
#' @param object object
#' @param N numbers of clones to for models
#' @param powerlaw_model_name residual of which powerlaw model to use?
#'   powerlaw_fixed/powerlaw_optim
#' @param snvs_name which snvs to to use?
#' @param method clustering method to use: BMix and mclust are currently supported.
#'   While mclust is a 3-4 times faster method, the BMix method is more accurate
#'   and usually fast enough.
#' @param upper_f_limit ignore variants with f higher than
#' @param verbose verbose?
#' @name fit_subclones
NULL



#' @rdname fit_subclones
#' @export
fit_subclones <- function(object,
                          N = 1:3,
                          powerlaw_model_name = active_models(object),
                          snvs_name = default_SNVs(object),
                          method = "BMix",
                          upper_f_limit = 0.75,
                          verbose = get_cevomod_verbosity()) {
  if (method == "BMix") {
    object <- object |>
      fit_subclones_bmix(N, powerlaw_model_name, snvs_name, upper_f_limit, verbose)
  } else if (method == "mclust") {
    object <- object |>
      fit_subclones_mclust(N, powerlaw_model_name, snvs_name, upper_f_limit, verbose)
  }

  object
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
  models <- get_models(object)
  expect_colnames <- c("N", "cellularity", "N_mutations")
  all(expect_colnames %in% colnames(models))
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
