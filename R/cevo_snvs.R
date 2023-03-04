
new_cevo_snvs <- function(tbl) {
  structure(tbl, class = c("cevo_snvs", class(tibble::tibble())))
}


#' Create cevo_snvs tibble
#' @param snvs tibble
#' @export
as_cevo_snvs <- function(snvs) {
  required_cols <- c("sample_id", "chrom", "pos", "VAF")
  optional_cols <- c("gene_symbol", "ref", "alt", "ref_reads", "alt_reads", "impact")

  missing_required_cols <- setdiff(required_cols, names(snvs))
  if (length(missing_required_cols)) {
    stop(str_c("input snvs is missing the following columns:", str_c(missing_required_cols, collapse = ", ")))
  }

  missing_optional_cols <- setdiff(optional_cols, names(snvs))
  optional_col_types <- list(
    gene_symbol = NA_character_, ref = NA_character_, alt = NA_character_,
    ref_reads = NA_real_, alt_reads = NA_real_, impact = NA_character_
  )
  for (col in missing_optional_cols) {
    snvs[[col]] <- optional_col_types[[col]]
  }

  snvs |>
    select(
      "sample_id", "chrom", "pos", "gene_symbol",
      "ref", "alt", "ref_reads", "alt_reads", "VAF", "impact",
      everything()
    ) |>
    new_cevo_snvs()
}



validate_SNVs <- function(snvs) {
  required_cols <- c(
    "sample_id", "chrom", "pos", "gene_symbol",
    "ref", "alt", "ref_reads", "alt_reads", "VAF", "impact"
  )
  missing_cols <- setdiff(required_cols, names(snvs))
  if (length(missing_cols)) {
    stop(str_c("snvs object is missing the following columns:", str_c(missing_cols, collapse = ", ")))
  }
}
