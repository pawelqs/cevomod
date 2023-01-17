
#' Get selected mutations based on 2 sample models
#' @param object object
#' @param ... other arguments
#' @export
identify_non_neutral_tail_mutations <- function(object, ...) {
  UseMethod("identify_non_neutral_tail_mutations")
}


#' @rdname identify_non_neutral_tail_mutations
#' @param sample1 sample1
#' @param sample2 sample2
#' @param method method
#' @param verbose lgl
#' @export
identify_non_neutral_tail_mutations.cevodata <- function(
        object,
        sample1 = NULL, sample2 = NULL,
        method = "basic",
        verbose = TRUE, ...) {
  if (!were_subclonal_models_fitted(object)) {
    stop("Fit subclonal models first!")
  }
  object[["joined_models"]] <- NULL

  splits <- object |>
    split_by("patient_id")
  splits <- splits |>
    map(
      identify_non_neutral_tail_mutations,
      sample1 = sample1, sample2 = sample2,
      method = method, verbose = verbose
    )
  object[["joined_models"]] <- splits |>
    map("joined_models") |>
    reduce(c)

  object
}


#' @rdname identify_non_neutral_tail_mutations
#' @export
identify_non_neutral_tail_mutations.singlepatient_cevodata <- function(
        object,
        sample1 = NULL, sample2 = NULL,
        method = "basic",
        verbose = TRUE, ...) {
  samples_data <- object$metadata |>
    select("sample", "sample_id")
  patient_id <- unique(object$metadata$patient_id)
  msg("Processing patient ", patient_id, "\t", new_line = FALSE, verbose = verbose)

  rowsample <- if (is.null(sample1)) samples_data$sample[[1]] else sample1
  colsample <- if (is.null(sample2)) samples_data$sample[[2]] else sample2
  sample_ids <- deframe(samples_data)[c(rowsample, colsample)]
  names(sample_ids) <- c("rows", "cols")

  mutations_mat <- object |>
    get_SNVs_2d_matrix(
      rows_sample = rowsample, cols_sample = colsample,
      verbose = verbose
    )
  zero_intervals <- c(rows = rownames(mutations_mat)[1], cols = colnames(mutations_mat)[1])
  row_predictions <- object |>
    get_binomial_model_predictions(
      sample_id = sample_ids["rows"], zero_interval = zero_intervals["rows"]
    )
  col_predictions <- object |>
    get_binomial_model_predictions(
      sample_id = sample_ids["cols"], zero_interval = zero_intervals["cols"]
    )

  if (method == "basic") {
    join_models <- solve_basic
  } else if (method == "MC") {
    join_models <- solve_MC
  }
  joined_models <- join_models(
    mutations_mat,
    row_predictions,
    col_predictions,
    N = 10, epochs = 1000, eps = 5,
    verbose = verbose
  )
  top_model_ev <- evaluate_MC_runs(joined_models$solution, row_predictions, col_predictions)

  object[["joined_models"]][[patient_id]] <- list(
    rowsample = rowsample,
    colsample = colsample,
    row_predictions = row_predictions,
    col_predictions = col_predictions,
    all_mutations_mat = mutations_mat,
    sel_mutations_mat = joined_models$solution,
    fit_metrics = top_model_ev,
    mc_arr = joined_models$mc_arr,
    metrics = joined_models$metrics,
    neutral_tail_mutations_mat = mutations_mat - joined_models$solution,
    selection_probability_mat = fill_na(joined_models$solution / mutations_mat, 0)
  )

  object[["models"]][["selection_probabilities"]] <- get_selection_probability_tbl(object)
  object
}


get_binomial_model_predictions <- function(object, sample_id, zero_interval = NULL) {
  binom_predictions <- get_residuals(object, model = "binomial_models") |>
    filter(.data$sample_id == .env$sample_id) |>
    select("VAF_interval", "binom_pred") |>
    deframe()
  zero_pred <- if (!is.null(zero_interval)) set_names(0, zero_interval) else NULL
  c(zero_pred, binom_predictions)
}


init_MC_simulation_limits <- function(mutations_mat, row_predictions, col_predictions) {
  upper_limits <- mutations_mat
  zero_interval <- rownames(mutations_mat)[[1]]
  for (row in rownames(mutations_mat)) {
    for (col in colnames(mutations_mat)) {
      if (row == zero_interval && col == zero_interval) {
        upper_limits[row, col] <- 0
      } else if (row == zero_interval) {
        upper_limits[row, col] <- min(mutations_mat[row, col], col_predictions[col] * 1.1)
      } else if (col == zero_interval) {
        upper_limits[row, col] <- min(mutations_mat[row, col], row_predictions[row] * 1.1)
      } else {
        upper_limits[row, col] <- min(
          mutations_mat[row, col],
          row_predictions[row] * 1.1,
          col_predictions[col] * 1.1
        )
      }
    }
  }
  upper_limits <- round(upper_limits)
  lower_limits <- upper_limits
  lower_limits[, ] <- 0
  limits <- list(lower = lower_limits, upper = upper_limits)
  class(limits) <- c("cevo_MC_sim_limits", class(limits))
  limits
}


run_MC_simulation <- function(upper_limits, lower_limits = NULL, iters = 2000) {
  mc_arr <- array(
    dim = c(nrow(upper_limits), ncol(upper_limits), iters),
    dimnames = list(rownames(upper_limits), colnames(upper_limits), 1:iters)
  )
  for (row in rownames(upper_limits)) {
    for (col in colnames(upper_limits)) {
      mc_arr[row, col, ] <- if (upper_limits[row, col] == 0) {
        0
      } else {
        stats::runif(
          n = iters,
          min = if (is.null(lower_limits)) 0 else lower_limits[row, col],
          max = upper_limits[row, col]
        )
      }
    }
  }
  mc_arr <- round(mc_arr)
  mc_arr
}


extract_models_to_tibble <- function(mc_arr, which) {
  which |>
    set_names(which) |>
    map(~mc_arr[, , .x]) |>
    map(as.data.frame) |>
    map(rownames_to_column, "rowsample") |>
    map(~pivot_longer(.x, -rowsample, names_to = "colsample", values_to = "N")) |>
    bind_rows(.id = "i")
}


average_solutions <- function(mc_arr, which = NULL) {
  mean_solution <- array(
    dim = dim(mc_arr)[1:2],
    dimnames = dimnames(mc_arr)[1:2]
  )
  if (is.null(which)) {
    which <- 1:dim(mc_arr)[[3]]
  }
  for (row in rownames(mean_solution)) {
    for (col in colnames(mean_solution)) {
      mean_solution[row, col] <- mean(mc_arr[row, col, which])
    }
  }
  mean_solution
  class(mean_solution) <- c("non_neutral_2d_fit", class(mean_solution))
  mean_solution
}


# ------------------------- method: basic ---------------------------------------

solve_basic <- function(mutations_mat, row_predictions, col_predictions,
                        N = 10, epochs = 1000, eps = 10,
                        verbose = TRUE) {
  limits <- init_MC_simulation_limits(mutations_mat, row_predictions, col_predictions)
  mc_arr <- run_MC_simulation(upper_limits = limits$upper, lower_limits = limits$lower, iters = N)

  for (i in 1:N) {
    msg(".", new_line = FALSE, verbose = verbose)
    mc_mat <- mc_arr[, , i]
    metrics <- evaluate_MC_runs(mc_mat, row_predictions, col_predictions)
    for (j in 1:epochs) {
      if (verbose == 2) {
        print(metrics$MSE)
      }

      ## New solution, randomize order of row and col fitting
      scaling_factors <- list(
          rows = calc_scaling_factors(rowSums(mc_mat), row_predictions)[-1],
          cols = calc_scaling_factors(colSums(mc_mat), col_predictions)[-1]
        ) |>
        map(enframe, name = "VAF_interval", value = "factor") |>
        bind_rows(.id = "dim") |>
        shuffle()

      for (k in seq_along(scaling_factors$dim)) {
        if (scaling_factors$dim[k] == "rows") {
          row <- scaling_factors$VAF_interval[k]
          new_values <- mc_mat[row, ] * scaling_factors$factor[k]
          mc_mat[row, ] <- pmin(new_values, limits$upper[row, ])
        } else if (scaling_factors$dim[k] == "cols") {
          col <- scaling_factors$VAF_interval[k]
          new_values <- mc_mat[, col] * scaling_factors$factor[k]
          mc_mat[, col] <- pmin(new_values, limits$upper[, col])
        }
      }

      ## Old solution
      # scaling_factors <- calc_scaling_factors(rowSums(mc_mat), row_predictions)
      # for (row in rownames(mc_mat)[-1]) {
      #   mc_mat[row, ] <- mc_mat[row, ] * scaling_factors[row]
      # }
      # mc_mat[mc_mat > limits$upper] <- limits$upper[mc_mat > limits$upper]
      #
      # scaling_factors <- calc_scaling_factors(colSums(mc_mat), col_predictions)
      # for (col in colnames(mc_mat)[-1]) {
      #   mc_mat[, col] <- mc_mat[, col] * scaling_factors[col]
      # }
      # mc_mat[mc_mat > limits$upper] <- limits$upper[mc_mat > limits$upper]

      mew_metrics <- evaluate_MC_runs(mc_mat, row_predictions, col_predictions)
      MSE_diff <- metrics$MSE - mew_metrics$MSE
      metrics <- mew_metrics
      mc_arr[, , i] <- mc_mat
      if (MSE_diff < eps) {
        msg("MSE_diff < eps, exiting", verbose = verbose == 2)
        break
      }
    }
  }
  msg("", verbose = verbose)

  metrics <- evaluate_MC_runs(mc_arr, row_predictions, col_predictions)
  final_solution <- average_solutions(mc_arr)

  list(
    solution = final_solution,
    mc_arr = mc_arr,
    metrics = metrics
  )
}


calc_scaling_factors <- function(actual, predicted) {
  scaling_factors <- predicted / actual
  scaling_factors[is.infinite(scaling_factors)] <- 1
  scaling_factors[is.na(scaling_factors)] <- 1
  scaling_factors
}


# ---------------------------- method: MC -------------------------------------

solve_MC <- function(mutations_mat, row_predictions, col_predictions, N = 5, verbose = TRUE) {
  limits <- init_MC_simulation_limits(mutations_mat, row_predictions, col_predictions)
  limit_ranges(limits) |> sum()
  for (i in 1:N) {
    limit_ranges(limits) |> sum() |> print()
    mc_arr <- run_MC_simulation(upper_limits = limits$upper, lower_limits = limits$lower, iters = 3000)
    metrics <- evaluate_MC_runs(mc_arr, row_predictions, col_predictions)
    # plot(metrics)

    top_models <- metrics |>
      filter(percent_rank(.data$MSE) < 0.01) |>
      arrange(percent_rank(.data$MSE))
    # plot_predictions_vs_fits(row_predictions, ev_res$rsums[top_models$i, ])
    # plot_predictions_vs_fits(col_predictions, ev_res$csums[top_models$i, ])

    mean_top_solution <- average_solutions(mc_arr, top_models$i)
    top_model_ev <- evaluate_MC_runs(mean_top_solution, row_predictions, col_predictions)
    top_model_ev |> print()

    limits <- tune_limits(limits, mc_arr, metrics, row_predictions, col_predictions, verbose = FALSE)
  }
  list(
    solution = mean_top_solution,
    mc_arr = mc_arr[, , top_models$i],
    metrics = metrics,
    top_models = top_models
  )
}


limit_ranges <- function(limits, zero.rm = TRUE) {
  zero_limits <- limits$upper == 0
  limits_diff <- limits$upper - limits$lower
  c(limits_diff[!zero_limits])
}


tune_limits <- function(limits, mc_arr, metrics, row_predictions, col_predictions, verbose = FALSE) {
  upper_limits <- limits$upper
  lower_limits <- limits$lower
  for (row in rownames(upper_limits)) {
    if (verbose) print(row)
    for (col in colnames(upper_limits)) {
      if (upper_limits[row, col] == lower_limits[row, col]) {
        next
      }
      metrics2 <- tibble(
        i = metrics$i,
        row_dist = (metrics$rsums[, row] - row_predictions[row])^2,
        col_dist = (metrics$csums[, col] - col_predictions[col])^2,
        dist = .data$row_dist + .data$col_dist
      )
      top_fits <- metrics2 |>
        mutate(rank = percent_rank(.data$dist)) |>
        filter(rank < 0.01)
      top_fits_arr <- mc_arr[, , top_fits$i]
      max_val <- max(top_fits_arr[row, col, ])
      min_val <- min(top_fits_arr[row, col, ])
      # mean_val <- mean(mc_arr[row, col, ])
      # sd <- sd(mc_arr[row, col, ])
      # max_val <- mean_val + 2 * sd
      # min_val <- mean_val - 2 * sd
      # upper
      if (max_val < row_predictions[row] || max_val < col_predictions[col]) {
        upper_limits[row, col] <- upper_limits[row, col]
      } else {
        upper_limits[row, col] <- max_val
      }
      # lower
      if (min_val > row_predictions[row] || min_val > col_predictions[col]) {
        # lower_limits[row, col] <- 0
        lower_limits[row, col] <- lower_limits[row, col]
      } else {
        lower_limits[row, col] <- min_val
      }
    }
  }
  upper_limits <- round(upper_limits)
  lower_limits <- round(lower_limits)
  limits <- list(lower = lower_limits, upper = upper_limits)
  class(limits) <- c("cevo_MC_sim_limits", class(limits))
  limits
}


# ----------------------- Models evaluation -----------------------------------

#' Evaluate Monte Carlo results
#' @param mc_arr Monte Carlo simutations array
#' @param rowsums_pred expected row sums
#' @param colsums_pred predicted col sums
#' @export
evaluate_MC_runs <- function(mc_arr, rowsums_pred, colsums_pred) {
  UseMethod("evaluate_MC_runs")
}


#' @export
evaluate_MC_runs.array <- function(mc_arr, rowsums_pred, colsums_pred) {
  iter <- dim(mc_arr)[[3]]
  rsums <- matrix(NA, nrow = iter, ncol = dim(mc_arr)[[1]])
  colnames(rsums) <- dimnames(mc_arr)[[1]]
  csums <- matrix(NA, nrow = iter, ncol = dim(mc_arr)[[2]])
  colnames(csums) <- dimnames(mc_arr)[[2]]
  rsums_MSE <- rep(NA_real_, iter)
  csums_MSE <- rep(NA_real_, iter)
  MSE <- rep(NA_real_, iter)
  # non_zero_sums <- matrix(NA, nrow = iter, ncol = length(row_predictions) - 1)
  for (i in 1:iter) {
    rsums[i, ] <- rowSums(mc_arr[, , i])
    csums[i, ] <- colSums(mc_arr[, , i])
    rsums_MSE[i] <- sum((rsums[i, -1] - rowsums_pred[-1])^2)
    csums_MSE[i] <- sum((csums[i, -1] - colsums_pred[-1])^2)
    MSE[i] <- rsums_MSE[i] + csums_MSE[i]
  }
  metrics <- tibble(
    i = 1:iter,
    rsums_MSE,
    csums_MSE,
    MSE,
    rsums,
    csums
  )
  class(metrics) <- c("cevo_MC_solutions_eval", class(metrics))
  metrics
}


#' @export
evaluate_MC_runs.matrix <- function(mc_arr, rowsums_pred, colsums_pred) {
  rsums <- rowSums(mc_arr)
  csums <- colSums(mc_arr)
  non_zero_sums <- c(rsums[-1], csums[-1])
  non_zero_preds <- c(rowsums_pred[-1], colsums_pred[-1])
  MSE <- sum((non_zero_sums - non_zero_preds)^2)
  rsq_rows <- rsq_vec(rowsums_pred[-1], rsums[-1])
  rsq_cols <- rsq_vec(colsums_pred[-1], csums[-1])
  rsq_total <- rsq_vec(non_zero_preds, non_zero_sums)
  res <- lst(MSE, rsq_rows, rsq_cols, rsq_total, rsums, csums)
  class(res) <- c("non_neutral_2d_fit_eval", class(res))
  res
}


#' @export
print.non_neutral_2d_fit_eval <- function(x, ...) {
  cli::cat_line("MSE:     ", round(x$MSE, digits = 2))
  cli::cat_line("Rows R^2:    ", round(x$rsq_rows, digits = 2))
  cli::cat_line("Cols R^2:    ", round(x$rsq_cols, digits = 2))
  cli::cat_line("Total R^2:   ", round(x$rsq_total, digits = 2))
  cli::cat_line("Row sums:  ", paste0(x$rsums[1:5], collapse = ", "), "...")
  cli::cat_line("Col sums:  ", paste0(x$csums[1:5], collapse = ", "), "...")
}


# @export
# print.cevo_MC_solutions_eval <- function(x, ...) {
#   MSE_range_str <- x$err$err |>
#     range() |>
#     round(digits = 0) |>
#     paste0(collapse = " - ")
#   cli::cat_line("Iters:     ", nrow(x$err))
#   cli::cat_line("MSE range: ", MSE_range_str)
# }


#' @export
plot.cevo_MC_solutions_eval <- function(x, ...) {
  graphics::hist(x$MSE)
}


# ----------------------- Get selection probabilities ------------------------


get_selection_probability_tbl <- function(object) {
  two_sample_probs <- get_single_sample_selection_probability(object)
  single_sample_probs <- get_two_sample_selection_probability(object)
  res <- two_sample_probs |>
    left_join(single_sample_probs, by = c("patient_id", "mutation_id"))
  res
}


get_two_sample_selection_probability <- function(object) {
  selection_prob_tbl <- object$joined_models |>
    map("selection_probability_mat") |>
    map(selection_prob_mat_to_long_tiblle) |>
    bind_rows(.id = "patient_id")
  samples <- colnames(selection_prob_tbl)[2:3]
  mutations <- get_SNVs_wider_intervals(object)

  res <- mutations |>
    left_join(selection_prob_tbl, by = c("patient_id", samples)) |>
    select("patient_id", "mutation_id", "2sample_selection_prob")
  res
}


selection_prob_mat_to_long_tiblle <- function(selection_probability_mat)  {
  rows_sample <- attributes(selection_probability_mat)$rows_sample
  cols_sample <- attributes(selection_probability_mat)$cols_sample
  selection_probability_mat |>
    as.data.frame() |>
    rownames_to_column(rows_sample) |>
    pivot_longer(-rows_sample, names_to = cols_sample, values_to = "2sample_selection_prob")
}


get_single_sample_selection_probability <- function(object) {
  binom_res <- get_residuals(object, model = "binomial_models") |>
    select("sample_id", "VAF_interval", "SFS", "binom_pred")
  metadata <- object$metadata[c("patient_id", "sample_id", "sample")]

  selection_prob_tbl <- metadata |>
    left_join(binom_res, by = "sample_id") |>
    transmute(
      .data$patient_id,
      .data$sample,
      .data$VAF_interval,
      `1sample_selection_prob` = binom_pred / SFS
    )

  SNVs(object) |>
    cut_VAF_intervals() |>
    unite_mutation_id() |>
    left_join(metadata, by = "sample_id") |>
    select("patient_id", "sample", "mutation_id", "VAF_interval") |>
    left_join(selection_prob_tbl, by = c("patient_id", "sample", "VAF_interval")) |>
    pivot_wider(names_from = "sample", values_from = "1sample_selection_prob") |>
    select("patient_id", "mutation_id", all_of(metadata$sample))
}


# -------------------------------- Plots --------------------------------------

plot_predictions_vs_fits <- function(predictions, fits) {
  plot(predictions)
  if (is.null(dim(fits))) {
    graphics::points(fits[-1], col = "red")
  } else {
    for (i in 1:nrow(fits)) {
      graphics::points(fits[i, -1], col = "red")
    }
  }
}

#' Plot heatmap of non neutral mutations
#' @param object object
#' @param ... other arguments
#' @export
plot_non_neutral_mutations_2D <- function(object, ...) {
  UseMethod("plot_non_neutral_mutations_2D")
}


#' @rdname plot_non_neutral_mutations_2D
#' @param colors vector of three colors to use
#' @export
plot_non_neutral_mutations_2D.cevodata <- function(object,
                                                   colors = c("black", "white", "red"),
                                                   ...) {
  joined_models <- object[["joined_models"]]
  joined_models |>
    map(function(x) {
      mat <- x$sel_mutations_mat
      plot_2d(
        mat,
        name = "N mutations",
        row_title = x$rowsample, column_title = x$colsample,
        col = circlize::colorRamp2(breaks = c(0, max(mat)/2, max(mat)), colors),
      )
    })
}


plot_2d <- function(mat, ...) {
  rlang::check_installed("ComplexHeatmap", reason = "Required to plot the heatmap")
  mat <- mat[rev(rownames(mat)), ]
  ComplexHeatmap::Heatmap(
    mat,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    ...
  )
}


#' @export
plot.non_neutral_2d_fit <- function(x, ...) {
  rlang::check_installed("pheatmap", reason = "to plot the 2d model")
  pheatmap::pheatmap(
    x, cluster_rows = FALSE, cluster_cols = FALSE
  )
}

