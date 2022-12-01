
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

  clones <- residuals |>
    nest_by(.data$sample_id) |>
    summarise(
      model = "binomial_clones",
      clones = fit_binomial_models(.data$data, N = N),
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
    discard(is.null) |>
    discard(any_clusters_overlap)

  if (length(mclust_res) > 0) {
    clones <- mclust_res |>
      map(mclust_to_clones_tbl, n_mutations = length(VAFs)) |>
      bind_rows() |>
      mutate(best = .data$BIC == max(.data$BIC))
  } else {
    clones <- empty_clones_tibble()
  }
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


empty_clones_tibble <- function() {
  tibble(
    n = integer(),
    component = character(),
    cellularity = double(),
    N_mutations = double(),
    BIC = double(),
    best = logical()
  )
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

