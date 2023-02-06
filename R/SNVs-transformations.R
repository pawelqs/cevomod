
#' Annotate mutation context and types for mutation signatures analysis
#' @param snvs snvs tbl
#' @param bsgenome BSGenome object
#' @export
annotate_mutation_contexts <- function(snvs, bsgenome) {
  rlang::check_installed("mutSignatures", reason = "Used to annotate mutations")
  snvs |>
    filter(str_length(.data$ref) == 1, str_length(.data$alt) == 1) |>
    mutSignatures::attachContext(
      chr_colName = "chrom",
      start_colName = "pos",
      end_colName = "pos",
      nucl_contextN = 3,
      BSGenomeDb = bsgenome
    ) |>
    mutSignatures::removeMismatchMut(
      refMut_colName = "ref",
      context_colName = "context",
      refMut_format = "N"
    ) |>
    mutSignatures::attachMutType(
      ref_colName = "ref",
      var_colName = "alt",
      context_colName = "context"
    )
}


#' Get matrix of per sample mutation types for mutation signatures analysis
#' @param snvs annotated snvs tbl
#' @return wide tbl
#' @export
count_mutation_types <- function(snvs) {
  rlang::check_installed("mutSignatures", reason = "Used to count mutations")
  counts <- mutSignatures::countMutTypes(
    snvs,
    mutType_colName = "mutType",
    sample_colName = "sample_id"
  )

  counts_mat <- counts@counts
  rownames(counts_mat) <- counts@mutTypes$mutTypes
  colnames(counts_mat) <- counts@sampleId$ID

  counts_mat |>
    rownames_to_column("MutationType")

}


get_SNVs_wider <- function(object, fill_na = NULL) {
  patients_to_samples <- object$metadata |>
    select("patient_id":"sample")

  snvs <- SNVs(object) |>
    select("sample_id", "chrom":"alt", "VAF") |>
    left_join(patients_to_samples, by = "sample_id") |>
    unite_mutation_id() |>
    select(-"sample_id") |>
    select("patient_id", everything()) |>
    pivot_wider(names_from = "sample", values_from = "VAF")

  if (!is.null(fill_na)) {
    snvs[is.na(snvs)] <- 0
  }
  snvs
}


unite_mutation_id <- function(snvs) {
  unite(snvs, "mutation_id", "chrom":"alt", sep = "-")
}


get_SNVs_2d_matrix <- function(object,
                               rows_sample = NULL, cols_sample = NULL,
                               bins = NULL, verbose = TRUE) {
  patients_to_samples <- object$metadata |>
    select("patient_id":"sample")

  if (is.null(rows_sample) || is.null(cols_sample)) {
    rows_sample <- patients_to_samples$sample[[1]]
    cols_sample <- patients_to_samples$sample[[2]]
    msg("Using '", rows_sample, "' as rows and '", cols_sample, "' as cols", verbose = verbose)
  }
  if (n_distinct(patients_to_samples$patient_id) > 1) {
    stop("This function works only for single sample objects")
  }
  if (rows_sample %not in% patients_to_samples$sample) {
    stop("Sample requested for rows is not present for this patient")
  }
  if (cols_sample %not in% patients_to_samples$sample) {
    stop("Sample requested for cols is not present for this patient")
  }

  breaks <- object |>
    SNVs() |>
    get_interval_breaks(bins = bins)
  sample_ids <- patients_to_samples$sample_id |>
    set_names(patients_to_samples$sample)
  rowsample_breaks <- breaks[[sample_ids[rows_sample]]]
  colsample_breaks <- breaks[[sample_ids[cols_sample]]]

  # Prepare SNVs wider tibble
  mutations <-get_SNVs_wider(object, fill_na = 0)
  mutations <- mutations[c("mutation_id", rows_sample, cols_sample)]
  colnames(mutations) <- c("mutation_id", "rows_sample", "cols_sample")

  # Cut intervals
  mutations <- mutations |>
    mutate(
      rows_sample = cut(rows_sample, breaks = rowsample_breaks),
      cols_sample = cut(cols_sample, breaks = colsample_breaks)
    )
  row_intervals <- levels(mutations$rows_sample)
  col_intervals <- levels(mutations$cols_sample)

  # Prepare square matrix
  incomplete_mat <- mutations |>
    mutate(across(c("rows_sample", "cols_sample"), as.character)) |>
    group_by(.data$rows_sample, .data$cols_sample) |>
    count() |>
    ungroup() |>
    pivot_wider(names_from = "cols_sample", values_from = "n") |>
    column_to_rownames("rows_sample")
  incomplete_mat[is.na(incomplete_mat)] <- 0

  mat <- matrix(0, nrow = length(row_intervals), ncol = length(col_intervals))
  rownames(mat) <- row_intervals
  colnames(mat) <- col_intervals
  mat[rownames(incomplete_mat), colnames(incomplete_mat)] <- as.matrix(incomplete_mat)

  attr(mat, "rows_sample") <- rows_sample
  attr(mat, "cols_sample") <- cols_sample
  attr(mat, "rows_sample_id") <- sample_ids[[rows_sample]]
  attr(mat, "cols_sample_id") <- sample_ids[[cols_sample]]
  mat
}


get_SNVs_wider_intervals <- function(object, fill_na = NULL, bins = NULL) {
  metadata <- object$metadata
  if (n_distinct(metadata$patient_id) > 1) {
    stop("This function works only for single sample objects")
  }

  breaks <- object |>
    SNVs() |>
    get_interval_breaks(bins = bins)
  breaks <- breaks[metadata$sample_id]
  names(breaks) <- metadata$sample

  mutations <-get_SNVs_wider(object, fill_na = 0)
  sample_intervals <- list()
  for (sample in metadata$sample) {
    VAF_intervals <- cut(mutations[[sample]], breaks = breaks[[sample]])
    sample_intervals[[sample]] <- levels(VAF_intervals)
    mutations[[sample]] <- as.character(VAF_intervals)
  }

  attr(mutations, "sample_intervals") <- sample_intervals
  mutations
}
