
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
  upper_limits <- round(upper_limits)[1:5, 1:5]
  # upper_limits[1:5, 1:5]

  iter <- 10000
  mc_arr <- array(
    dim = c(nrow(mutations_mat), ncol(mutations_mat), iter),
    dimnames = list(rownames(mutations_mat), colnames(mutations_mat), 1:iter)
  )
  for (row in rownames(mutations_mat)) {
    for (col in colnames(mutations_mat)) {
      mc_arr[row, col, ] <- round(runif(iter, max = upper_limits[row, col]))
    }
  }
  mc_arr[, , 1]

  # Evaluate
  err <- rep(NA_real_, iter)
  rsums <- matrix(NA, nrow = iter, ncol = nrow(mutations_mat))
  csums <- matrix(NA, nrow = iter, ncol = ncol(mutations_mat))
  # non_zero_sums <- matrix(NA, nrow = iter, ncol = length(row_predictions) - 1)
  for (i in 1:iter) {
    rsums[i, ] <- rowSums(mc_arr[, , i])
    csums[i, ] <- colSums(mc_arr[, , i])
    non_zero_sums <- c(rsums[i, -1], csums[i, -1])
    non_zero_preds <- c(row_predictions[-1], col_predictions[-1])
    err[i] <- sum((non_zero_sums - non_zero_preds)^2)
  }
  hist(err)

  err <- tibble(err, i = 1:iter) |>
    mutate(rank = percent_rank(err))

  selected_solutions <- err |>
    filter(rank < 0.01) |>
    arrange(err)

  sel <- err$i[[1]]

  plot(row_predictions)
  points(rsums[sel, -1], col = "red")


  plot(col_predictions)
  points(csums[sel, -1], col = "red")

  pheatmap::pheatmap(
    mc_arr[, , sel], cluster_rows = FALSE, cluster_cols = FALSE
  )
}



