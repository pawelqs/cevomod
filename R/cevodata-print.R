
#' @export
print.cevodata <- function(x, ...) {
  summ <- summary(x)
  SNV_assays_str <- if (length(summ$SNV_assays) > 0) {
    summ$SNV_assays[summ$SNV_assays == summ$active_SNVs] <- paste0(summ$active_SNVs, " (default)")
    summ$SNV_assays |>
      str_c(collapse = ", ")
  } else {
    "None"
  }
  CNV_assays_str <- if (length(summ$CNV_assays) > 0) {
    summ$CNV_assays[summ$CNV_assays == summ$active_CNVs] <- paste0(summ$active_CNVs, " (default)")
    summ$CNV_assays |>
      str_c(collapse = ", ")
  } else {
    "None"
  }

  cli::cat_line("<cevodata> dataset: ", x$name, col = "#45681e")
  cli::cat_line("Genome: ", summ$genome, col = "#628f2f")
  cli::cat_line("SNV assays: ", SNV_assays_str)
  cli::cat_line("CNV assays: ", CNV_assays_str)
  cli::cat_line(summ$metadata_str)
  cli::cat_line(summ$SNVs_str)
  cli::cat_line("Active models: ", x$active_models)
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
      select("patient_id", "sample_id") |>
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
  snvs_added <- object$active_SNVs != ""
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
      select("patient_id", "chrom":"alt") |>
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
