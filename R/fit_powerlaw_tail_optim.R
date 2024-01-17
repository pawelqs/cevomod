# -------------------------- fit_powerlaw_tail_optim --------------------------

#' Fitting power-law tails with aptimum exponent value
#'
#' `fit_powerlaw_tail_optim()` uses `stats::optim` to find optimal A and alpha
#' whch maximizes SFS area under the powerlaw curve (*sampled* region of SFS and
#' the range of f values below the maximum SFS value does not count) and
#' minimizes negative error - where the curve is above the real SFS (*sampled*
#' are does not count). Penalty for the negative error depends on the number of
#' points with the negative error value. Penalty value is the sum of error values
#' to the power of x, where x is the length of the vector of negative error values.
#' This allows the powerlaw curve to detach from the SFS for 1 or two bins, but
#' then the penalty rises extremely.
#'
#' @param object cevodata
#' @param name name in the models' slot
#' @param av_filter average filter values to be applied to f
#' @param control control param of stats::optim()
#' @param verbose verbose?
#' @param ... other arguments passed to stats::optim()
#' @examples
#' data("tcga_brca_test")
#' cd <- tcga_brca_test |>
#'   dplyr::filter(sample_id %in% c("TCGA-AC-A23H-01","TCGA-AN-A046-01")) |>
#'   fit_powerlaw_tail_optim()
#' @name powerlaw_optim
NULL


#' @rdname powerlaw_optim
#' @export
fit_powerlaw_tail_optim <- function(object, ...) {
  UseMethod("fit_powerlaw_tail_optim")
}


#' @rdname powerlaw_optim
#' @inheritParams get_non_zero_SFS_range
#' @param peak_detection_upper_limit Upper f value up to which the main peak is searched
#' @param reward_upper_limit Mutations under the curve up to this limit will be rewarded
#' @param bootstraps Number of bootstrap samples, or FALSE to make no resampling.
#'   **This option significantly extendis the model fitting time!!**
#' @export
fit_powerlaw_tail_optim.cevodata <- function(object,
                                             name = "powerlaw_optim",
                                             bootstraps = FALSE,
                                             allowed_zero_bins = 2,
                                             y_treshold = 1,
                                             y_threshold_pct = 0.01,
                                             av_filter = c(1/3, 1/3, 1/3),
                                             peak_detection_upper_limit = 0.3,
                                             reward_upper_limit = 0.4,
                                             control = list(maxit = 1000, ndeps = c(0.1, 0.01)),
                                             verbose = get_verbosity(),
                                             ...) {
  if (!bootstraps) {
    sfs <- get_SFS(object, name = "SFS")
    coefs <- fit_powerlaw_tail_optim(
      sfs,
      name = name,
      allowed_zero_bins = allowed_zero_bins,
      y_treshold = y_treshold,
      y_threshold_pct = y_threshold_pct,
      av_filter = av_filter,
      peak_detection_upper_limit = peak_detection_upper_limit,
      reward_upper_limit = reward_upper_limit,
      control = control,
      verbose = verbose
    )
    residuals <- coefs |>
      filter(.data$best) |>
      calc_powerlaw_model_residuals(sfs)
    info <- list(f_column = attr(sfs, "f_column"))

    models <- lst(coefs, residuals, info)
    class(models) <- c("cevo_powerlaw_models", "list")
    add_models(object, models, name = name)
    # object$models[[name]] <- models
    # object <- calc_powerlaw_model_residuals(object, name)
    # object$active_models <- name
    # object
  } else {
    rlang::check_installed("rsample", reason = "to perform bootstrap sampling of SNVs")
    sfs_resamples <- calc_SFS_resamples(object, times = bootstraps, verbose = verbose)

    models <- sfs_resamples |>
      map(
        fit_powerlaw_tail_optim,
        name = name,
        allowed_zero_bins = allowed_zero_bins,
        y_treshold = y_treshold,
        y_threshold_pct = y_threshold_pct,
        av_filter = av_filter,
        peak_detection_upper_limit = peak_detection_upper_limit,
        reward_upper_limit = reward_upper_limit,
        control = control,
        verbose = verbose
    )

    bootstrap_models <- models |>
      map("bootstrap_models") |>
      bind_rows()

    models <- models |>
      map("models") |>
      bind_rows(.id = "sample_id")

    bootstrap_name <- str_c(name, "_bootstraps")
    object$models[[bootstrap_name]] <- bootstrap_models
    object <- calc_powerlaw_model_residuals(object, bootstrap_name)
    object$models[[name]] <- models
    object <- calc_powerlaw_model_residuals(object, name)
    object$active_models <- name
    # object
    add_models(object, models, name = name)
  }
}


#' @rdname powerlaw_optim
#' @export
fit_powerlaw_tail_optim.cevo_SFS_bootstraps <- function(object,
                                                        name = "powerlaw_optim",
                                                        allowed_zero_bins = 2,
                                                        y_treshold = 1,
                                                        y_threshold_pct = 0.01,
                                                        av_filter = c(1/3, 1/3, 1/3),
                                                        peak_detection_upper_limit = 0.3,
                                                        reward_upper_limit = 0.4,
                                                        control = list(maxit = 1000, ndeps = c(0.1, 0.01)),
                                                        verbose = get_verbosity(),
                                                        ...) {
  rlang::check_installed("rsample", reason = "to perform bootstrap sampling of SNVs")
  msg("Fitting models to ", unique(object$sfs[[1]]$sample_id), " resamples", verbose = verbose)
  pb <- if (verbose) progress_bar$new(total = nrow(object)) else NULL

  object$models <- object$sfs |>
    map(function(sfs) {
      if (!is.null(pb)) pb$tick()
      fit_powerlaw_tail_optim(
        sfs,
        name = name,
        allowed_zero_bins = allowed_zero_bins,
        y_treshold = y_treshold,
        y_threshold_pct = y_threshold_pct,
        av_filter = av_filter,
        peak_detection_upper_limit = peak_detection_upper_limit,
        reward_upper_limit = reward_upper_limit,
        control = control,
        verbose = verbose_down(verbose)
      )
    })

  object$tidy_models <- object$models |>
    map(~pivot_longer(.x, all_of(c("A", "alpha")), names_to = "term", values_to = "estimate"))

  conf_intervals <- rsample::int_pctl(object, .data$tidy_models)

  bootstrap_models <- object$models |>
    set_names(object$id) |>
    bind_rows(.id = "resample_id")
  class(bootstrap_models) <- c("cevo_bootstrap_powerlaw_models", class(bootstrap_models))

  models <- conf_intervals |>
    pivot_wider(
      names_from = "term",
      values_from = ".lower":".upper",
      names_glue = "{term}{.value}",
      names_vary = "slowest"
    ) |>
    transmute(
      model = name,
      component = "powerlaw tail",
      A = .data$A.estimate,
      .data$A.lower, .data$A.upper,
      alpha = .data$alpha.estimate,
      .data$alpha.lower, .data$alpha.upper
    )
  class(models) <- c("cevo_powerlaw_models", class(models))

  lst(bootstrap_models, models)
}


#' @rdname powerlaw_optim
#' @export
fit_powerlaw_tail_optim.cevo_SFS_tbl <- function(object,
                                                 name = "powerlaw_optim",
                                                 allowed_zero_bins = 2,
                                                 y_treshold = 1,
                                                 y_threshold_pct = 0.01,
                                                 av_filter = c(1/3, 1/3, 1/3),
                                                 peak_detection_upper_limit = 0.3,
                                                 reward_upper_limit = 0.4,
                                                 control = list(maxit = 1000, ndeps = c(0.1, 0.01)),
                                                 verbose = get_verbosity(),
                                                 ...) {
  msg("Fitting optimized power-law models...", verbose = verbose)
  start_time <- Sys.time()

  sfs <- object
  bounds <- sfs |>
    get_non_zero_SFS_range(
      allowed_zero_bins = allowed_zero_bins,
      y_treshold = y_treshold,
      y_threshold_pct = y_threshold_pct
    ) |>
    rename(lower_bound = "from", higher_bound = "to")
  # bounds <- get_f_range(SNVs(object), pct_left = pct_left, pct_right = pct_right)
  nbins <- summarise(sfs, nbins = n(), .by = "sample_id")

  data <- sfs |>
    left_join(bounds, by = "sample_id") |>
    filter(.data$f > .data$lower_bound, .data$f < .data$higher_bound) |>
    mutate(length = n(), .by = "sample_id") |>
    filter(length >= 3) |>
    select("sample_id", "f", "y") |>
    mutate(
      y = .data$y |>
        stats::filter(av_filter) |>
        fill_na(0),
      .by = "sample_id"
    ) |>
    nest_by(.data$sample_id) |>
    left_join(nbins, by = "sample_id") |>
    expand_grid(
      init_A = c(1, 2, 4, 8, 16, 32),
      init_alpha = c(0.8, 1.2, 1.8, 2.5, 3.5)
    )

  pb <- if (verbose) progress_bar$new(total = nrow(data)) else NULL
  models <- data |>
    rowwise("sample_id") |>
    summarise(
      model = name,
      component = "powerlaw tail",
      opt = td_optim(
        .data$init_A, .data$init_alpha,
        .data$data,
        peak_detection_upper_limit = peak_detection_upper_limit,
        reward_upper_limit = reward_upper_limit,
        control = control,
        pb = pb,
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
  models
}


# --------------------------------- Helpers ------------------------------------

td_optim <- function(init_A, init_alpha, data,
                     peak_detection_upper_limit = 0.3,
                     reward_upper_limit = 0.4,
                     control, pb = NULL, ...) {
  res <- stats::optim(
    par = c(init_A, init_alpha),
    fn = td_objective_function,
    x = data$f,
    y = data$y,
    peak_detection_upper_limit = peak_detection_upper_limit,
    reward_upper_limit = reward_upper_limit,
    control = control,
    ...
  )

  if (!is.null(pb)) pb$tick()
  res
}


td_objective_function <- function(params, x, y,
                                  peak_detection_upper_limit = 0.3,
                                  reward_upper_limit = 0.4
                                  ) {
  A <- params[[1]]
  alpha <- params[[2]]
  y1 <- A * 1/(x ^ alpha)
  err <- y - y1
  err[err < -1000000] <- -1000000

  err_segments <- segment(err < 0)
  sampled_range <- err < 0 & err_segments == 0

  # Reward for number of mutations under the curve
  y2 <- pmin(y1, y)
  before_max <- seq_along(y) < which.max(y[x < peak_detection_upper_limit]) # detects peak up to f = 0.3
  y2[before_max | sampled_range | (x > reward_upper_limit)] <- 0
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
    mutate(best = row_number() == 1) |>
    ungroup()
}


#' Get range of non empty SFS bins
#' @param sfs SFS
#' @param allowed_zero_bins number of allowed empty bins in the interval
#' @param y_treshold bins with less mutations will be considered empty
#' @param y_threshold_pct bins that have less mutations than this param times the
#'   height of the higherst peak will be considered empty
#' @keywords internal
get_non_zero_SFS_range <- function(sfs,
                                   allowed_zero_bins = 1,
                                   y_treshold = 1,
                                   y_threshold_pct = 0.01) {
  sfs |>
    group_by(.data$sample_id) |>
    mutate(
      empty_bin = (.data$y < y_treshold) | (.data$y < max(.data$y) * y_threshold_pct),
      segment_number = segment(.data$empty_bin),
      empty_low_f_range = .data$empty_bin & .data$segment_number == 0
    ) |>
    filter(!.data$empty_low_f_range) |>
    group_by(.data$sample_id, .data$segment_number) |>
    mutate(
      segment_length = n(),
      keep = !.data$empty_bin | (.data$empty_bin & (.data$segment_length <= allowed_zero_bins))
    ) |>
    group_by(.data$sample_id) |>
    mutate(new_segments = segment(.data$keep)) |>
    filter(.data$new_segments == 0) |>
    summarise(
      from = min(.data$f),
      to = max(.data$f)
    )
}
