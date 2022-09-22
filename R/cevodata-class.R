
new_cevodata <- function(name, genome) {
  cd <- list(
    name = name,
    genome = genome,
    SNVs = list(),
    CNVs = list(),
    clones = list(),
    models = list(),
    active_SNVs = NULL,
    active_CNVs = NULL,
    active_clones = NULL
  )
  structure(cd, class = "cevodata")
}


#' Create new cevomod dataset object
#'
#' @param name dataset name
#' @param genome genome version
#' @param snvs tibble with SNVs
#' @param snvs_name name for SNVs assay
#' @param cnvs tibble with CNVs
#' @param cnvs_name name for CNVs assay
#' @return `cevodata` object
#'
#' @export
init_cevodata <- function(name, genome = NULL,
                          snvs = NULL, snvs_name = NULL,
                          cnvs = NULL, cnvs_name = NULL) {
  genome <- if (is.null(genome)) "unknown" else genome
  cd <- new_cevodata(name, genome)
  if (!is.null(snvs)) {
    cd <- add_SNV_data(cd, snvs, snvs_name)
  }
  if (!is.null(cnvs)) {
    cd <- add_CNV_data(cd, cnvs, cnvs_name)
  }
  cd
}


#' @export
print.cevodata <- function(x, ...) {
  summ <- summary(x)
  SNV_assays_str <- if (length(summ$SNV_assays) > 0) {
    summ$SNV_assays |>
      str_c(collapse = ", ") |>
      str_replace(summ$active_SNVs, str_c(summ$active_SNVs, " (default)"))
  } else {
    "None"
  }
  CNV_assays_str <- if (length(summ$CNV_assays) > 0) {
    summ$CNV_assays |>
      str_c(collapse = ", ") |>
      str_replace(summ$active_CNVs, str_c(summ$active_CNVs, " (default)"))
  } else {
    "None"
  }

  cli::cat_line("<cevodata> dataset: ", x$name, col = "#45681e")
  cli::cat_line("Genome: ", summ$genome, col = "#628f2f")
  cli::cat_line("SNV assays: ", SNV_assays_str)
  cli::cat_line("CNV assays: ", CNV_assays_str)
  cli::cat_line(summ$n_patients, " patients, ", summ$n_samples, " samples")
  if (!is.null(x$active_SNVs)) {
    cli::cat_line("SNVs:")
    print(x$SNVs[[x$active_SNVs]])
  }
  if (!is.null(x$active_CNVs)) {
    cli::cat_line("CNVs:")
    print(x$CNVs[[x$active_CNVs]])
  }
}


#' @export
summary.cevodata <- function(object, ...) {
  list(
    name = object$name,
    genome = object$genome,
    SNV_assays = names(object$SNVs),
    active_SNVs = object$active_SNVs,
    CNV_assays = names(object$CNVs),
    active_CNVs = object$active_CNVs,
    n_patients = n_distinct(object$SNVs[[object$active_SNVs]]$patient_id),
    n_samples = n_distinct(object$SNVs[[object$active_SNVs]]$sample_id)
  )
}


#' @return tibble
#' @describeIn assays Get SNVs from cevodata dataset
#' @export
SNVs.cevodata <- function(object, which = object$active_SNVs, ...) {
  if (which %not in% names(object$SNVs)) {
    stop(str_c(which, " does not exist in object$SNVs"))
  }
  object$SNVs[[which]]
}


#' @describeIn assays Add new SNVs to cevodata
#' @export
add_SNV_data.cevodata <- function(object, snvs, name = NULL, ...) {
  if(is.null(name)) {
    n <- length(object$SNVs)
    name <- if (n == 0) "snvs" else str_c("snvs", n)
  }
  validate_SNVs(snvs)
  object$SNVs[[name]] <- snvs
  default_SNVs(object) <- name
  object
}


validate_SNVs <- function(snvs) {
  required_cols <- c(
    "patient_id", "sample_id", "sample", "chrom", "pos", "gene_symbol",
    "ref", "alt", "ref_reads", "alt_reads", "VAF", "impact"
  )
  missing_cols <- setdiff(required_cols, names(snvs))
  if (length(missing_cols)) {
    stop(str_c("snvs object is missing the following columns:", str_c(missing_cols, collapse = ", ")))
  }
}


#' @describeIn active_assays Get default SNVs assay of cevodata
#' @export
default_SNVs.cevodata <- function(object, ...) {
  object$active_SNVs
}


#' @describeIn active_assays Set default SNVs assay of cevodata
#' @export
`default_SNVs<-.cevodata` <- function(object, ..., value) {
  if (value %not in% names(object$SNVs)) {
    stop("Chosen SNV assay must exist in object$SNVs")
  }
  object$active_SNVs <- value
  object
}


#' @describeIn assays Get CNVs from cevodata dataset
#' @export
CNVs.cevodata <- function(object, which = object$active_CNVs, ...) {
  if (which %not in% names(object$CNVs)) {
    stop(str_c(which, " does not exist in object$CNVs"))
  }
  object$CNVs[[which]]
}


#' @describeIn assays Add new CNVs to cevodata
#' @export
add_CNV_data.cevodata <- function(object, cnvs, name = NULL, ...) {
  if(is.null(name)) {
    n <- length(object$CNVs)
    name <- if (n == 0) "cnvs" else str_c("cnvs", n)
  }
  validate_CNVs(cnvs)
  object$CNVs[[name]] <- cnvs
  default_CNVs(object) <- name
  object
}


validate_CNVs <- function(cnvs) {
  required_cols <- c(
    "patient_id", "sample_id", "sample", "chrom", "start", "end",
    "log_ratio", "BAF", "total_cn", "major_cn", "minor_cn"
  )
  missing_cols <- setdiff(required_cols, names(cnvs))
  if (length(missing_cols)) {
    stop(str_c("cnvs object is missing the following columns:", str_c(missing_cols, collapse = ", ")))
  }
}


#' @describeIn active_assays Get default CNVs assay of cevodata
#' @export
default_CNVs.cevodata <- function(object, ...) {
  object$active_CNVs
}


#' @describeIn active_assays Set default CNVs assay of cevodata
#' @export
`default_CNVs<-.cevodata` <- function(object, ..., value) {
  if (value %not in% names(object$CNVs)) {
    stop("Chosen CNV assay must exist in object$CNVs")
  }
  object$active_CNVs <- value
  object
}
