
#' readthis integration
#'
#' @description
#' [readthis](https://github.com/pawelqs/readthis) package may be used to easily
#' read the data from some popular mutation callers into R environment. readthis
#' functions can be supplied not only with the single file paths, but also with
#' lists of files or even paths to the directories with files to be loaded (and
#' cevodata object is to store the data from many samples!)
#'
#' readthis functions return tibbles or list of tibbles. These tibbles/
#' objects usually are instances of *cevo_<software_name>* S3 classes. cevomod
#' implements methods that allow to add these types of data to the cevodata
#' objects conveniently.
#'
#' @param cd <cevodata> object
#' @param data Object read with readthis functions
#' @param name Name for the data
#' @param verbose Verbose?
#' @param ... Other arguments
#'
#' @examples
#' # library(cevomod)
#'
#' ascat_dir <- system.file("extdata", "ASCAT", package = "readthis")
#' ascat <- readthis::read_ascat_files(ascat_dir)
#' cd <- init_cevodata("Test dataset") |>
#'   add_data(ascat)
#'
#' @name readthis-integration
NULL



#' @describeIn readthis-integration add_data() function takes cevodata as the
#'   first argument, so it is a preferred method for adding data in R pipelines.
#' @export
add_data <- function(cd, data, ...) {
  add_to_cevodata(data, cd)
}


#' @describeIn readthis-integration add_to_cevodata() is a generic with a set
#'   of methods for different classes of `data`. These methods are called by
#'   add_data() function.
#' @export
add_to_cevodata <- function(data, cd, name, verbose, ...) {
  UseMethod("add_to_cevodata")
}


#' @export
add_to_cevodata.cevo_ASCAT <- function(data, cd,
                                       name = "ASCAT",
                                       verbose = get_cevomod_verbosity(),
                                       ...) {
  sample_data <- data$sample_statistics |>
    mutate(ascat_purity = 1 - .data$normal_contamination)
  cd |>
    add_CNV_data(data$cnvs, name = name) |>
    add_sample_data(sample_data) |>
    use_purity("ascat_purity", verbose = verbose)
}


#' @export
add_to_cevodata.cevo_FACETS <- function(data, cd,
                                        name = "FACETS",
                                        verbose = get_cevomod_verbosity(),
                                        ...) {
  cnvs <- data |>
    select(-"Purity", -"Ploidy")
  sample_data <- data |>
    select("sample_id", facets_purity = "Purity", facets_ploidy = "Ploidy") |>
    unique()
  cd |>
    add_CNV_data(data, name = name) |>
    add_sample_data(sample_data) |>
    use_purity("facets_purity", verbose = verbose)
}


#' @export
add_to_cevodata.cevo_Mutect <- function(data, cd,
                                        name = "Mutect",
                                        verbose = get_cevomod_verbosity(),
                                        ...) {
  patient_ids_present <- "patient_id" %in% names(data)

  if (patient_ids_present) {
    sample_data <- data |>
      select("patient_id", "sample_id") |>
      unique()
    data$patient_id <- NULL
  }

  cd <- add_SNV_data(cd, data, name = name)
  if (patient_ids_present) {
    cd <- add_sample_data(cd, sample_data)
  }

  cd
}


#' @export
add_to_cevodata.cevo_Strelka <- function(data, cd,
                                         name = "Strelka",
                                         verbose = get_cevomod_verbosity(),
                                         ...) {
  patient_ids_present <- "patient_id" %in% names(data)

  if (patient_ids_present) {
    sample_data <- data |>
      select("patient_id", "sample_id") |>
      unique()
    data$patient_id <- NULL
  }

  cd <- add_SNV_data(cd, data, name = name)
  if (patient_ids_present) {
    cd <- add_sample_data(cd, sample_data)
  }

  cd
}
