
new_cevodata <- function(name, genome) {
  cd <- list(
    name = name,
    genome = genome,
    metadata = NULL,
    SNVs = list(),
    CNVs = list(),
    clones = list(),
    models = list(),
    active_SNVs = "",
    active_CNVs = "",
    active_clones = "",
    active_model = ""
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


#' Get/Add SNV/CNV data from the cevodata dataset
#' @param object object
#' @param snvs tibble with SNVs
#' @param cnvs tibble with CNVs
#' @param name name for SNVs/CNVs assay
#' @param which assay to use - uses active_SNVs if nonei
#' @param ... other arguments
#' @name assays


#' @rdname assays
#' @export
SNVs <- function(object, ...) {
  UseMethod("SNVs")
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


#' @rdname assays
#' @export
add_SNV_data <- function(object, ...) {
  UseMethod("add_SNV_data")
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


#' @rdname assays
#' @export
CNVs <- function(object, ...) {
  UseMethod("CNVs")
}

#' @describeIn assays Get CNVs from cevodata dataset
#' @export
CNVs.cevodata <- function(object, which = object$active_CNVs, ...) {
  if (which %not in% names(object$CNVs)) {
    stop(str_c(which, " does not exist in object$CNVs"))
  }
  object$CNVs[[which]]
}


#' @rdname assays
#' @export
add_CNV_data <- function(object, ...) {
  UseMethod("add_CNV_data")
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


#' Get/Set active assays of the cevodata object
#' @param object object
#' @param value name of new default assay
#' @param ... other arguments
#' @name active_assays


#' @rdname active_assays
#' @export
default_SNVs <- function(object, ...) {
  UseMethod("default_SNVs")
}

#' @describeIn active_assays Get default SNVs assay of cevodata
#' @export
default_SNVs.cevodata <- function(object, ...) {
  object$active_SNVs
}


#' @rdname active_assays
#' @export
`default_SNVs<-` <- function(object, ..., value) {
  UseMethod("default_SNVs<-")
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


#' @rdname active_assays
#' @export
default_CNVs <- function(object, ...) {
  UseMethod("default_CNVs")
}

#' @describeIn active_assays Get default CNVs assay of cevodata
#' @export
default_CNVs.cevodata <- function(object, ...) {
  object$active_CNVs
}

#' @rdname active_assays
#' @export
`default_CNVs<-` <- function(object, ..., value) {
  UseMethod("default_CNVs<-")
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


#' Get names of CNV variables
#' @param object object
#' @param ... other arguments
#' @export
get_CNVs_var_names <- function(object, ...) {
  UseMethod("get_CNVs_var_names")
}

#' @describeIn get_CNVs_var_names Get CNV variable names from cevodata object
#' @param which CNV assay to use
#' @export
get_CNVs_var_names.cevodata <- function(object, which = default_CNVs(object), ...) {
  cnvs_metadata <- CNVs(object, which = which) |>
    select(-.data$sample_id, -.data$chrom, -.data$start, -.data$end) |>
    drop_na_columns()
  colnames(cnvs_metadata)
}


#' Add metadata to the cevodata object
#' @param object object
#' @param data name of new default assay
#' @param ... other arguments
#' @name cevo_metadata


#' @rdname cevo_metadata
#' @export
add_patient_data <- function(object, ...) {
  UseMethod("add_patient_data")
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


#' @rdname cevo_metadata
#' @export
add_sample_data <- function(object, ...) {
  UseMethod("add_sample_data")
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
