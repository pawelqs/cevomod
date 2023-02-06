
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
  bounds <- get_VAF_range(SNVs(object), pct_left = 0.02, pct_right = 0.98)
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
      opt = optim(
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

  zeroed_count <- 0
  for (i in seq_along(err)) {
    if (err[[i]] < 0) {
      err[[i]] <- 0
      zeroed_count <- zeroed_count + 1
    }
    if (err[[i]] > 0)
      break
  }

  tail_muts <- pmin(y, y1)[-(1:zeroed_count)]
  # tail_mut_count <- sum(tail_muts) / sum(y)
  weights <- (1 - x[-(1:zeroed_count)]) ^ 2
  # weights <- (1 - x) ^ 2
  # weights <- (length(tail_muts):1) / length(tail_muts)^2
  mut_reward <- sum(tail_muts*weights)
  # mut_reward <- sum((tail_muts*weights)^2)

  scaled_err <- err
  scaled_err[is.infinite(scaled_err)] <- 0
  negatives <- scaled_err < 0
  segmented <- segment(negatives)
  segments <- unique(segmented)

  weights <- segmented
  for (seg in segments) {
    n <- sum(segmented == seg)
    weights[segmented == seg] <- n
  }
  scaled_err[scaled_err > 0] <- 0
  scaled_err <- scaled_err / y
  scaled_err[is.infinite(scaled_err)] <- 0
  scaled_err <- abs(scaled_err) + 1
  # too_high_penalty <- sum(scaled_err[scaled_err < 0] ^ 2, na.rm = TRUE)
  too_high_penalty <- sum(scaled_err^weights)

  maxim <- mut_reward - too_high_penalty
  -maxim
}
