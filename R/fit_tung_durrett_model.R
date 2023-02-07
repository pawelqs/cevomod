
#' Fitting Tung Durrett models
#'
#' @param object cevodata
#' @param name name in the models' slot
#' @param verbose verbose?
#' @param ... other arguments
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
                                             verbose = TRUE, ...) {
  msg("Fitting Tung-Durrett models...", verbose = verbose)
  sfs <- get_SFS(object, name = "SFS")
  bounds <- get_VAF_range(SNVs(object), pct_left = 0, pct_right = 0.98)
  nbins <- get_sample_sequencing_depths(SNVs(object)) |>
    transmute(.data$sample_id, nbins = .data$median_DP)

  data <- sfs |>
    left_join(bounds, by = "sample_id") |>
    filter(.data$VAF > .data$lower_bound, .data$VAF < .data$higher_bound) |>
    select("sample_id", "VAF", "y") |>
    nest_by(.data$sample_id) |>
    left_join(nbins, by = "sample_id")

  models <- data |>
    summarise(
      model = "tung_durrett",
      component = "powerlaw",
      opt = stats::optim(
        par = c(5, 1),
        fn = td_objective_function,
        x = data$VAF,
        y = data$y
      ) |> list(),
      A = .data$opt$par[[1]] * round(.data$nbins),
      alpha = .data$opt$par[[2]],
      val = -.data$opt$value,
      .groups = "drop"
    ) |>
    select(-"opt")
  class(models) <- c("cevo_powerlaw_models", class(models))

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
  y2[before_max | sampled_range] <- 0
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
