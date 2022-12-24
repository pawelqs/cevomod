
# Get selected mutations based on 2 sample models
get_selected_mutations <- function(object, ...) {
  # As for now for one patient only
  if (!were_subclonal_models_fitted(object)) {
    stop("Fit subclonal models first!")
  }

  samples_data <- object$metadata |>
    select(sample_id:sample)
  rowsample <- samples_data$sample[[1]]
  colsample <- samples_data$sample[[2]]

  mutations_mat <- object |>
    get_SNVs_2d_matrix(rows_sample = rowsample, cols_sample = colsample, bins = 100)

  intervals <- tibble(
    interval = rownames(mutations_mat),
    centers = get_interval_centers(interval)
  )

  VAF_zero_pred <- tibble(
    VAF = 0, binom_pred = 0,
    sample = samples_data$sample
  )
  binom_predictions <- get_residuals(object) |>
    select(sample_id, VAF, binom_pred) |>
    left_join(samples_data, by = "sample_id") |>
    select(-sample_id)
  binom_predictions <- VAF_zero_pred |>
    bind_rows(binom_predictions) |>
    pivot_wider(names_from = "sample", values_from = "binom_pred") |>
    select(all_of(c("VAF", rowsample, colsample))) |>
    set_names(c("VAF", "rowsample", "colsample"))

  predictions_by_interval <- binom_predictions |>
    rebinarize_distribution(VAFs = intervals$centers) |>
    bind_cols(intervals)
  row_predictions <- deframe(predictions_by_interval[, c("interval", "rowsample")])
  col_predictions <- deframe(predictions_by_interval[, c("interval", "colsample")])

  join_models <- if (method == "basic") {
    solve_basic
  } else if (method == "MC") {
    solve_MC
  }
  x <- join_models(mutations_mat, row_predictions, col_predictions, N = 10, epochs = 1000, eps = 5)
  # Loop
  plot_predictions_vs_fits(row_predictions, x$metrics$rsums)
  plot_predictions_vs_fits(col_predictions, x$metrics$csums)

  plot_predictions_vs_fits(row_predictions, rowSums(x$solution))
  plot_predictions_vs_fits(col_predictions, colSums(x$solution))
  plot(x$solution)
  top_model_ev <- evaluate_MC_runs(x$solution, row_predictions, col_predictions)
  top_model_ev
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
        upper_limits[row, col] <- min(mutations_mat[row, col], row_predictions[row] * 1.1, col_predictions[col] * 1.1)
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
        runif(
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

solve_basic <- function(mutations_mat, row_predictions, col_predictions, N = 10, epochs = 1000, eps = 10) {
  limits <- init_MC_simulation_limits(mutations_mat, row_predictions, col_predictions)

  mc_arr <- run_MC_simulation(upper_limits = limits$upper, lower_limits = limits$lower, iters = N)

  for (i in 1:N) {
    mc_mat <- mc_arr[, , i]
    metrics <- evaluate_MC_runs(mc_mat, row_predictions, col_predictions)
    for (j in 1:epochs) {
      print(metrics$MSE)

      rsums <- rowSums(mc_mat)
      scaling_factors <- row_predictions / rsums
      scaling_factors[is.infinite(scaling_factors)] <- 1
      scaling_factors[is.na(scaling_factors)] <- 1
      for (row in rownames(mc_mat)[-1]) {
        mc_mat[row, ] <- mc_mat[row, ] * scaling_factors[row]
      }
      mc_mat[mc_mat > limits$upper] <- limits$upper[mc_mat > limits$upper]

      csums <- colSums(mc_mat)
      scaling_factors <- col_predictions / csums
      scaling_factors[is.infinite(scaling_factors)] <- 1
      scaling_factors[is.na(scaling_factors)] <- 1
      for (col in colnames(mc_mat)[-1]) {
        mc_mat[, col] <- mc_mat[, col] * scaling_factors[col]
      }
      mc_mat[mc_mat > limits$upper] <- limits$upper[mc_mat > limits$upper]

      mew_metrics <- evaluate_MC_runs(mc_mat, row_predictions, col_predictions)
      MSE_diff <- metrics$MSE - mew_metrics$MSE
      metrics <- mew_metrics
      mc_arr[, , i] <- mc_mat
      if (MSE_diff < eps) {
        message("MSE_diff < eps, exiting")
        break
      }
    }
  }
  metrics <- evaluate_MC_runs(mc_arr, row_predictions, col_predictions)
  final_solution <- average_solutions(mc_arr)
  list(
    solution = final_solution,
    mc_arr = mc_arr,
    metrics = metrics
  )
}


# ------------------------- method: MC ---------------------------------------

solve_MC <- function(mutations_mat, row_predictions, col_predictions, N = 5) {
  limits <- init_MC_simulation_limits(mutations_mat, row_predictions, col_predictions)
  limit_ranges(limits) |> sum()
  for (i in 1:N) {
    limit_ranges(limits) |> sum() |> print()
    mc_arr <- run_MC_simulation(upper_limits = limits$upper, lower_limits = limits$lower, iters = 3000)
    metrics <- evaluate_MC_runs(mc_arr, row_predictions, col_predictions)
    # plot(metrics)

    top_models <- metrics |>
      filter(percent_rank(MSE) < 0.01) |>
      arrange(percent_rank(MSE))
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
        dist = row_dist + col_dist
      )
      top_fits <- metrics2 |>
        mutate(rank = percent_rank(dist)) |>
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
  hist(x$MSE)
}


# -------------------------------- Plots --------------------------------------

plot_predictions_vs_fits <- function(predictions, fits) {
  plot(predictions)
  if (is.null(dim(fits))) {
    points(fits[-1], col = "red")
  } else {
    for (i in 1:nrow(fits)) {
      points(fits[i, -1], col = "red")
    }
  }
}


#' @export
plot.non_neutral_2d_fit <- function(x, ...) {
  pheatmap::pheatmap(
    x, cluster_rows = FALSE, cluster_cols = FALSE
  )
}

