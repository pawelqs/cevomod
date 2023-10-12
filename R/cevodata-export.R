

## ------------------------------ CliP --------------------------------

#' Export cevodata to CliP input
#'
#' [CliP](https://github.com/wwylab/CliP) is an algorithm for clonal structure
#' identification through penalizing pairwise differences by Wenyi Wang Lab
#' at MD Anderson Cancer Center in Houston.
#'
#' @param cd cevodata object
#' @param out_dir Directory where to save files. List of items is returned if
#'   out_dir is NULL
#' @param snvs_name name of the snvs to use
#' @param cnvs_name name of the cnvs to use
#' @param purity_column name of the metadata column with the purity estimates
#'   to be used
#' @param keep_chromosomes list of non-sex chromosomes. CliP does not use sex
#'   chromosomes
#' @export
to_clip <- function(cd, out_dir = NULL,
                    snvs_name = default_SNVs(cd),
                    cnvs_name = default_CNVs(cd),
                    purity_column = "purity",
                    keep_chromosomes = str_c("chr", 1:22)) {
  snvs <- get_clip_snvs(cd, snvs_name, keep_chromosomes)
  cnvs <- get_clip_cnvs(cd, cnvs_name, keep_chromosomes)
  purities <- get_clip_purities(cd, purity_column)

  clip_data <- lst(snvs, cnvs, purities) |>
    transpose()

  if (is.null(out_dir)) {
    clip_data
  } else {
    save_clip_files(clip_data, out_dir)
  }
}



get_clip_snvs <- function(cd,
                          snvs_name = default_SNVs(cd),
                          keep_chromosomes = str_c("chr", 1:22)) {
  empty_clip_snvs <- tibble(
    chromosome_index = double(),
    position = integer(),
    alt_count = integer(),
    ref_count = integer()
  )

  SNVs(cd, which = snvs_name) |>
    filter(.data$chrom %in% keep_chromosomes) |>
    transmute(
      sample_id = parse_factor(.data$sample_id, levels = get_sample_ids(cd)),
      chromosome_index = chromosomes_to_int(.data$chrom),
      position = .data$pos,
      alt_count = .data$alt_reads,
      ref_count = .data$ref_reads
    ) |>
    nest_by(.data$sample_id) |>
    complete(.data$sample_id, fill = list(data = list(empty_clip_snvs))) |>
    deframe()
}



get_clip_cnvs <- function(cd,
                          cnvs_name = default_CNVs(cd),
                          keep_chromosomes = str_c("chr", 1:22)) {
  empty_clip_cnvs <- tibble(
    chromosome_index = double(),
    start_position = double(),
    end_position = double(),
    major_cn = double(),
    minor_cn = double(),
    total_cn = double()
  )

  CNVs(cd, which = cnvs_name) |>
    filter(.data$chrom %in% keep_chromosomes) |>
    transmute(
      sample_id = parse_factor(.data$sample_id, levels = get_sample_ids(cd)),
      chromosome_index = chromosomes_to_int(.data$chrom),
      start_position = .data$start,
      end_position = .data$end,
      major_cn = .data$major_cn,
      minor_cn = .data$minor_cn,
      total_cn = .data$total_cn
    ) |>
    filter(.data$total_cn > 0) |>   # these records break CliP
    nest_by(.data$sample_id) |>
    complete(.data$sample_id, fill = list(data = list(empty_clip_cnvs))) |>
    deframe()
}



get_clip_purities <- function(cd, purity_column = "purity") {
  cd$metadata |>
    select("sample_id", all_of(purity_column)) |>
    deframe()
}



chromosomes_to_int <- function(chrom) {
  case_when(
    chrom %in% str_c("chr", 1:22) ~ str_replace(chrom, "chr", "") |> as.integer(),
    chrom == "chrX" ~ 23,
    chrom == "chrY" ~ 24,
    chrom == "chrMT" ~ 25,
    TRUE ~ NA_integer_
  )
}



save_clip_files <- function(clip_data, out_dir) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }

  files <- imap(
    clip_data,
    function(x, sample_id) {
      # sample_id2 <- str_replace(sample_id, " ", "_")
      snvs_file <- file.path(out_dir, str_c(sample_id, ".snv.tsv"))
      cnvs_file <- file.path(out_dir, str_c(sample_id, ".cnv.tsv"))
      purity_file <- file.path(out_dir, str_c(sample_id, ".purity.tsv"))

      write_tsv(x$snvs, snvs_file)
      write_tsv(x$cnvs, cnvs_file)
      write_file(as.character(x$purities), purity_file)

      tibble(sample_id, snvs_file, cnvs_file, purity_file)
    }
  )

  invisible(bind_rows(files))
}


