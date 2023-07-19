
new_cevodata <- function(name, genome, cancer) {
  cd <- list(
    name = name,
    genome = genome,
    cancer = cancer,
    metadata = NULL,
    SNVs = list(),
    CNVs = list(),
    models = list(),
    misc = list(),
    misc_by_patient = list(),
    misc_by_sample = list(),
    active_SNVs = "",
    active_CNVs = "",
    active_models = ""
  )
  structure(cd, class = "cevodata")
}


#' Create new cevomod dataset object
#'
#' @param name dataset name
#' @param genome genome version
#' @param cancer cancer type from `driver_genes` tbl
#' @param snvs tibble with SNVs
#' @param snvs_name name for SNVs assay
#' @param cnvs tibble with CNVs
#' @param cnvs_name name for CNVs assay
#' @return `cevodata` object
#'
#' @export
init_cevodata <- function(name, genome = "unknown", cancer = "unknown",
                          snvs = NULL, snvs_name = NULL,
                          cnvs = NULL, cnvs_name = NULL) {
  cd <- new_cevodata(name, genome, cancer)
  if (!is.null(snvs)) {
    cd <- add_SNV_data(cd, snvs, snvs_name)
  }
  if (!is.null(cnvs)) {
    cd <- add_CNV_data(cd, cnvs, cnvs_name)
  }
  cd
}


#' Set cancer type for the object
#' @param object object to set cancer type
#' @param ... other arguments
#' @export
set_cancer_type <- function(object, ...) {
  UseMethod("set_cancer_type")
}

#' @describeIn set_cancer_type Set cancer type for cevodata object
#' @param cancer_type cancer type
#' @export
set_cancer_type.cevodata <- function(object, cancer_type, ...) {
  object$cancer <- cancer_type
  object
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
  snvs <- as_cevo_snvs(snvs)
  object$SNVs[[name]] <- snvs
  default_SNVs(object) <- name
  meta <- snvs |>
    select("sample_id") |>
    unique() |>
    as_tibble()
  object <- add_sample_data(object, meta)
  object
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
    select("sample_id") |>
    unique() |>
    as_tibble()
  object <- add_sample_data(object, meta)
  object
}


validate_CNVs <- function(cnvs) {
  required_cols <- c(
    "sample_id", "chrom", "start", "end"
    # "log_ratio", "BAF", "total_cn", "major_cn", "minor_cn"
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
      select("patient_id", "sample_id", "sample", everything())
  }
  object
}


#' Choose purity measure
#'
#' <cevodata> metadata can contain purity measures in columns other than 'purity'.
#' T his function can be used to set 'purity' values using values from requested
#' column
#'
#' @param cd <cevodata> object
#' @param name Name of the metadata column with chosen purity values
#' @param verbose Verbose?
#' @export
use_purity <- function(cd, name, verbose = get_cevomod_verbosity()) {
  if (name %not in% names(cd$metadata)) {
    stop(
      "`name` should be a name of the column in the metadata tibble, ",
      "which should be used as purity measure"
    )
  } else {
    msg("Using '", name, "' as default purity measure", verbose = verbose)
    if (!is.null(cd$metadata[["purity"]])) {
      cd$metadata$prev_purity <- cd$metadata$purity
    }
    cd$metadata$purity <- cd$metadata[[name]]
    cd
  }
}


is_cevodata_singlepatient <- function(object) {
  n_patients <- count_patients(object)
  if (is.na(n_patients)) {
    return(FALSE)
  } else if (n_patients > 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


count_patients <- function(cd) {
  if (!is.null(cd$metadata[["patient_id"]])) {
    n_distinct(cd$metadata$patient_id)
  } else {
    NA_integer_
  }
}
