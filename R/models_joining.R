
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
    get_2d_SNVs_matrix(rows_sample = rowsample, cols_sample = colsample, bins = 100)

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

  iter <- 2000
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




}


get_2d_SNVs_matrix <- function(object,
                               rows_sample = NULL, cols_sample = NULL,
                               bins = 100) {
  patients_to_samples <- object$metadata |>
    select(patient_id:sample)

  if (is.null(rows_sample) || is.null(cols_sample)) {
    rows_sample <- patients_to_samples$sample[[1]]
    cols_sample <- patients_to_samples$sample[[2]]
    message("Using '", rows_sample, "' as rows and '", cols_sample, "' as cols")
  }
  if (rows_sample %not in% patients_to_samples$sample) {
    stop("Sample requested for rows is not present for this sample")
  }
  if (cols_sample %not in% patients_to_samples$sample) {
    stop("Sample requested for cols is not present for this sample")
  }

  breaks <- c(-0.1, seq(0, 1, length.out = bins + 1))

  # Prepare SNVs wider tibble
  mutations <-get_SNVs_wider(object, fill_na = 0)
  mutations <- mutations[c("mutation_id", rows_sample, cols_sample)]
  colnames(mutations) <- c("mutation_id", "rows_sample", "cols_sample")

  # Cut intervals
  mutations <- mutations |>
    mutate(across(where(is.double), ~cut(.x, breaks = breaks)))
  intervals <- levels(mutations$rows_sample)

  # Prepare square matrix
  incomplete_mat <- mutations |>
    mutate(across(c("rows_sample", "cols_sample"), as.character)) |>
    group_by(rows_sample, cols_sample) |>
    count() |>
    ungroup() |>
    pivot_wider(names_from = cols_sample, values_from = n) |>
    column_to_rownames("rows_sample")
  incomplete_mat[is.na(incomplete_mat)] <- 0

  mat <- matrix(0, nrow = length(intervals), ncol = length(intervals))
  colnames(mat) <- intervals
  rownames(mat) <- intervals
  mat[rownames(incomplete_mat), colnames(incomplete_mat)] <- as.matrix(incomplete_mat)

  attr(mat, "rows_sample") <- rows_sample
  attr(mat, "cols_sample") <- cols_sample
  mat
}


get_interval_centers <- function(intervals) {
  intervals |>
    str_replace_all("[\\(\\)\\[\\]]", "") |>
    str_split(pattern = ",") |>
    map(parse_double) |>
    map(mean) |>
    unlist()
}


get_interval_width <- function(intervals) {
  intervals |>
    str_replace_all("[\\(\\)\\[\\]]", "") |>
    str_split(pattern = ",") |>
    map(parse_double) |>
    map(~.x[[2]] - .x[[1]]) |>
    unlist() |>
    median()
}

