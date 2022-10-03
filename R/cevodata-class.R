
new_cevodata <- function(name, genome) {
  cd <- list(
    name = name,
    genome = genome,
    metadata = NULL,
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
  cli::cat_line(summ$metadata_str)
  cli::cat_line(summ$SNVs_str)
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
  summ <- list(
    name = object$name,
    genome = object$genome,
    SNV_assays = names(object$SNVs),
    active_SNVs = object$active_SNVs,
    CNV_assays = names(object$CNVs),
    active_CNVs = object$active_CNVs
  )
  summ <- c(
    summ,
    summarize_metadata(object),
    summarize_SNVs(object)
  )
  summ
}


summarize_metadata <- function(object) {
  meta <- object$metadata
  patient_id_present <- !is.null(meta[["patient_id"]])

  samples_per_patient <- if (patient_id_present) {
    meta |>
      select(.data$patient_id, .data$sample_id) |>
      unique() |>
      group_by(.data$patient_id) |>
      count()
  }

  summ <- list(
    n_patients = if (patient_id_present) n_distinct(meta$patient_id),
    n_samples = n_distinct(meta$sample_id),
    samples_per_patient_min = if (patient_id_present) min(samples_per_patient$n),
    samples_per_patient_max = if (patient_id_present) max(samples_per_patient$n)
  )
  summ[["metadata_str"]] <- if (patient_id_present) {
    samples_per_patient <- stringify_number_of_samples(
      summ$samples_per_patient_min,
      summ$samples_per_patient_max
    )
    str_c(summ$n_patients, " cases, ", summ$n_samples, " samples, ", samples_per_patient)
  } else {
    str_c(summ$n_samples, " samples")
  }
  summ
}


stringify_number_of_samples <- function(min, max) {
  if (is.null(min) && is.null(max)) {
    "No samples added"
  } else if (min == 1 && max == 1) {
    "1 sample per case"
  } else if (min == max) {
    str_c(min, " samples per case")
  } else {
    str_c(min, " - ", max, " samples per case")
  }
}


summarize_SNVs <- function(object) {
  snvs_added <- !is.null(object$active_SNVs)
  patient_id_present <- !is.null(object$metadata[["patient_id"]])

  summ <- list()
  if (snvs_added) {
    mutations_per_sample <- SNVs(object) |>
      group_by(.data$sample_id) |>
      count()
    summ$mutations_per_sample_mean = round(mean(mutations_per_sample$n))
    summ$mutations_per_sample_sd = round(sd(mutations_per_sample$n))
    summ$n_mutations <- nrow(SNVs(object))
    summ$SNVs_str <- str_c(
      summ$n_mutations, " mutations total, ",
      summ$mutations_per_sample_mean, " +/- ", summ$mutations_per_sample_sd,
      " mutations per sample"
    )
  }
  if (patient_id_present) {
    mutations_per_patient <- SNVs(object) |>
      left_join(object$metadata, by = "sample_id") |>
      select(.data$patient_id, .data$chrom:.data$alt) |>
      unique() |>
      group_by(.data$patient_id) |>
      count()
    summ$mutations_per_patient_mean = round(mean(mutations_per_patient$n))
    summ$mutations_per_patient_sd = round(sd(mutations_per_patient$n))
    summ$SNVs_str <- str_c(
      summ$n_mutations, " mutations total, ",
      summ$mutations_per_patient_mean, " +/- ", summ$mutations_per_patient_sd,
      " mutations per case"
    )
  }
  if (!snvs_added & !patient_id_present) {
    summ$SNVs_str <- str_c("No mutations added")
  }
  summ
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
  meta <- snvs |>
    select(.data$sample_id) |>
    unique() |>
    as_tibble()
  object <- add_sample_data(object, meta)
  object
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
  meta <- cnvs |>
    select(.data$sample_id) |>
    unique() |>
    as_tibble()
  object <- add_sample_data(object, meta)
  object
}


validate_CNVs <- function(cnvs) {
  required_cols <- c(
    "sample_id", "chrom", "start", "end",
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


#' @describeIn cevo_metadata Add patient data to cevodata object
#' @export
add_patient_data.cevodata <- function(object, data, ...) {
  if ("patient_id" %not in% colnames(data)) {
    stop("Data must have patient_id column!")
  }
  if (is.null(object$metadata)) {
    object$metadata <- data
  } else {
    keys <- intersect(c("patient_id", "sample_id", "sample"), colnames(data))
    object$metadata <- full_join(object$metadata, data, by = keys)
  }
  object
}


#' @describeIn cevo_metadata Add samples' data to cevodata object
#' @export
add_sample_data.cevodata <- function(object, data, ...) {
  if ("sample_id" %not in% colnames(data)) {
    stop("Data must have sample_id column!")
  }
  if (is.null(object$metadata)) {
    object$metadata <- data
  } else {
    keys <- c("patient_id", "sample_id", "sample") |>
      intersect(colnames(data)) |>
      intersect(colnames(object$metadata))
    object$metadata <- full_join(object$metadata, data, by = keys)
  }
  if (all(c("patient_id", "sample_id", "sample") %in% colnames(object$metadata))) {
    object$metadata <- object$metadata |>
      select(.data$patient_id, .data$sample_id, .data$sample, everything())
  }
  object
}


#' Filter/subset cevodata object
#'
#' This is a wrapper around dplyr::filter function which can be used to subset
#' cevodata object. Works like dplyr::filter, performs the filtering on metadata,
#' then filters SNVs, CNVs, clones and models keeping samples kept in metadata
#'
#' @param .data cevodata object
#' @param ... expression passed to dplyr::filter(...)
#' @inheritParams dplyr::filter
#' @return cevodata object
#' @export
filter.cevodata <- function(.data, ..., .preserve = FALSE) {
  new_object <- .data
  new_object$metadata <- filter(new_object$metadata, ...)
  ids <- new_object$metadata$sample_id
  new_object$SNVs <- map(new_object$SNVs, ~filter(.x, sample_id %in% ids))
  new_object$CNVs <- map(new_object$CNVs, ~filter(.x, sample_id %in% ids))
  new_object$clones <- map(new_object$clones, ~filter(.x, sample_id %in% ids))
  new_object$models <- map(new_object$models, ~filter(.x, sample_id %in% ids))
  new_object
}


#' Merge two cevodata objects
#' @inheritParams base::merge
#' @param name Name of the merged object
#' @export
merge.cevodata <- function(x, y, name = "Merged datasets", ...) {
  genome <- if (x$genome == y$genome) x$genome else "multiple genomes"
  metadata <- bind_rows(x$metadata, y$metadata)
  cd <- init_cevodata(name = name, genome = genome) |>
    add_sample_data(metadata)

  cd$SNVs <- bind_assays(x, y, "SNVs")
  cd$CNVs <- bind_assays(x, y, "CNVs")
  cd$clones <- bind_assays(x, y, "clones")
  cd$models <- bind_assays(x, y, "models")

  message("Setting active SNVs to ", x$active_SNV)
  message("Setting active CNVs to ", x$active_CNV)
  message("Setting active clones to ", x$active_clones)
  cd$active_SNVs <- NULL
  cd$active_CNVs <- NULL
  cd$active_clones <- NULL
  active_assays <- list(
    active_SNVs = x$active_SNVs,
    active_CNVs = x$active_CNVs,
    active_clones = x$active_clones
  )
  cd <- c(cd, active_assays)
  structure(cd, class = "cevodata")
}


bind_assays <- function(x, y, slot_name) {
  assays <- c(names(x[[slot_name]]), names(y[[slot_name]])) |>
    unique()
  if (length(assays) > 0) {
    assays |>
      set_names(assays) |>
      map(~bind_rows(x[[slot_name]][[.x]], y[[slot_name]][[.x]]))
  } else {
    list()
  }
}
