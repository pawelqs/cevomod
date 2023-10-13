
## --------------------------- cevo_snvs class --------------------------------

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
      "ref", "alt", "ref_reads", "alt_reads", "impact", "VAF",
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


## ---------------------------- Misc functions --------------------------------


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


#' Unite many columns to create mutation_id column
#' @param snvs SNVs
#' @param sep Separator
#' @param remove Remove united columns?
#' @export
unite_mutation_id <- function(snvs, sep = "-", remove = TRUE) {
  unite(snvs, "mutation_id", "chrom":"alt", sep = sep, remove = remove)
}


#' Filter SNVs by position: using regions tbl or bed file
#' @param snvs snvs tbl with columns: sample_id, chrom, pos
#' @param regions regions tbl with columns chrom, start, end
#' @param bed_file bed file
#' @export
filter_SNVs_by_regions <- function(snvs, regions = NULL, bed_file = NULL) {
  if (is.null(regions) && is.null(bed_file)) {
    stop("Provide one of: regions, bed_file")
  }

  if (!is.null(bed_file)) {
    # bed_file = "tests/testdata/regions.tsv" # for tests only
    regions <- bed_file |>
      read_tsv(col_types = "cii", col_names = c("chrom", "start", "end"))

    # Unlike the coordinate system used by other standards such as GFF, the system
    # used by the BED format is zero-based for the coordinate start and one-based
    # for the coordinate end.
    regions <- regions |>
      mutate(start = .data$start + 1)
  }
  snv_classes <- class(snvs)

  regions_gr <- regions |>
    rename(seqnames = "chrom") |>
    plyranges::as_granges()
  snvs_gr <- snvs |>
    mutate(
      seqnames = .data$chrom,
      start = .data$pos,
      end = .data$pos,
      .before = "sample_id"
    ) |>
    plyranges::as_granges()

  filtered_snvs <- plyranges::filter_by_overlaps(snvs_gr, regions_gr) |>
    as_tibble() |>
    select(-("seqnames":"strand"))
  class(filtered_snvs) <- snv_classes

  filtered_snvs
}
