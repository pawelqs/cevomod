
#' Fitting Tung Durrett models
#'
#' `fit_tung_durrett_models()` uses `stats::optim` to find optimal A and alpha
#' whch maximizes SFS area under the powerlaw curve (*sampled* region of SFS and
#' the range of VAF values below the maximum SFS value does not count) and
#' minimizes negative error - where the curve is above the real SFS (*sampled*
#' are does not count). Penalty for the negative error depends on the number of
#' points with the negative error value. Penalty value is the sum of error values
#' to the power of x, where x is the length of the vector of negative error values.
#' This allows the powerlaw curve to detach from the SFS for 1 or two bins, but
#' then the penalty rises extremely.
#'
#' @param object cevodata
#' @param name name in the models' slot
#' @param pct_left drop pct of the lowest frequency mutations
#' @param pct_right drop pct of the highest frequency mutations
#' @param control control param of stats::optim()
#' @param verbose verbose?
#' @param ... other arguments passed to stats::optim()
#' @examples
#' data("tcga_brca_test")
#' cd <- tcga_brca_test |>
#'   dplyr::filter(sample_id %in% c("TCGA-AC-A23H-01","TCGA-AN-A046-01")) |>
#'   fit_tung_durrett_models()
#' @name tung_durrett
NULL


#' @rdname tung_durrett
#' @export
fit_tung_durrett_models <- function(object, ...) {
  UseMethod("fit_tung_durrett_models")
}


#' @rdname tung_durrett
#' @export
fit_tung_durrett_models.cevodata <- function(object,
                                             name = "tung_durrett",
                                             pct_left = 0, pct_right = 0.98,
                                             control = list(maxit = 1000, ndeps = c(0.1, 0.01)),
                                             verbose = TRUE, ...) {
  msg("Fitting Tung-Durrett models...", verbose = verbose)
  start_time <- Sys.time()

  sfs <- get_SFS(object, name = "SFS")
  # bounds <- get_non_zero_SFS_range(sfs, y_treshold = 2, allowed_zero_bins = 2) |>
  #   rename(lower_bound = "from", higher_bound = "to")
  bounds <- get_VAF_range(SNVs(object), pct_left = pct_left, pct_right = pct_right)
  nbins <- get_sample_sequencing_depths(SNVs(object)) |>
    transmute(.data$sample_id, nbins = .data$median_DP)

  data <- sfs |>
    left_join(bounds, by = "sample_id") |>
    filter(.data$VAF > .data$lower_bound, .data$VAF < .data$higher_bound) |>
    select("sample_id", "VAF", "y") |>
    nest_by(.data$sample_id) |>
    left_join(nbins, by = "sample_id") |>
    expand_grid(
      init_A = c(1, 2, 4, 8, 16, 32),
      init_alpha = c(0.8, 1.2, 1.8, 2.5, 3.5)
    ) |>
    rowwise("sample_id")

  models <- data |>
    summarise(
      model = "tung_durrett",
      component = "powerlaw",
      opt = stats::optim(
        par = c(.data$init_A, .data$init_alpha),
        fn = td_objective_function,
        x = data$VAF,
        y = data$y,
        control = control,
        ...
      ) |> list(),
      A = .data$opt$par[[1]] * round(.data$nbins),
      alpha = .data$opt$par[[2]],
      convergence = .data$opt$convergence,
      message = .data$opt$message,
      value = -.data$opt$value,
      .groups = "drop"
    ) |>
    select(-"opt") |>
    evaluate_td_models() |>
    filter(.data$best)
  class(models) <- c("cevo_powerlaw_models", class(models))

  msg("Models fitted in ", Sys.time() - start_time, " seconds", verbose = verbose)

  object$models[[name]] <- models
  object <- calc_powerlaw_model_residuals(object, "tung_durrett")
  object$active_models <- name
  object
}


td_objective_function <- function(params, x, y) {
  A <- params[[1]]
  alpha <- params[[2]]
  y1 <- A * 1/(x ^ alpha)
  err <- y - y1
  err[err < -1000000] <- -1000000

  err_segments <- segment(err < 0)
  sampled_range <- err < 0 & err_segments == 0

  # Reward for number of mutations under the curve
  y2 <- pmin(y1, y)
  before_max <- seq_along(y) < which.max(y[x < 0.3]) # detects peak up to VAF = 0.3
  y2[before_max | sampled_range | (x > 0.4)] <- 0
  mut_reward <- sum(y2)

  # Penalty for bins too low for the curve
  weights <- err_segments
  for (seg in err_segments) {
    seg_length <- sum(err_segments == seg)
    weights[err_segments == seg] <- seg_length
  }

  too_low_err <- err
  too_low_err[sampled_range] <- 0
  too_low_err[too_low_err > 0] <- 0
  too_low_err <- abs(too_low_err)
  too_high_penalty <- sum(too_low_err^weights)

  maxim <- mut_reward - too_high_penalty
  -maxim
}


evaluate_td_models <- function(tbl) {
  tbl |>
    group_by(.data$sample_id) |>
    arrange(desc(.data$value), .by_group = TRUE) |>
    mutate(best = row_number() == 1)
}
