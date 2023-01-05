
get_SNVs_wider <- function(object, fill_na = NULL) {
  patients_to_samples <- object$metadata |>
    select("patient_id":"sample")

  snvs <- SNVs(object) |>
    select("sample_id", "chrom":"alt", "VAF") |>
    left_join(patients_to_samples, by = "sample_id") |>
    unite("mutation_id", "chrom":"alt", sep = "-") |>
    select(-"sample_id") |>
    select("patient_id", everything()) |>
    pivot_wider(names_from = "sample", values_from = "VAF")

  if (!is.null(fill_na)) {
    snvs[is.na(snvs)] <- 0
  }
  snvs
}


get_SNVs_2d_matrix <- function(object,
                               rows_sample = NULL, cols_sample = NULL,
                               bins = 100) {
  patients_to_samples <- object$metadata |>
    select("patient_id":"sample")

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
    group_by(.data$rows_sample, .data$cols_sample) |>
    count() |>
    ungroup() |>
    pivot_wider(names_from = "cols_sample", values_from = "n") |>
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
