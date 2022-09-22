
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


#' @rdname assays
#' @export
add_SNV_data <- function(object, ...) {
  UseMethod("add_SNV_data")
}

#' @rdname assays
#' @export
CNVs <- function(object, ...) {
  UseMethod("CNVs")
}


#' @rdname assays
#' @export
add_CNV_data <- function(object, ...) {
  UseMethod("add_CNV_data")
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


#' @rdname active_assays
#' @export
`default_SNVs<-` <- function(object, ..., value) {
  UseMethod("default_SNVs<-")
}


#' @rdname active_assays
#' @export
default_CNVs <- function(object, ...) {
  UseMethod("default_CNVs")
}


#' @rdname active_assays
#' @export
`default_CNVs<-` <- function(object, ..., value) {
  UseMethod("default_CNVs<-")
}


#' Site Frequency Spectra
#'
#' Creates  cevodata$models$SFS with the groupping variables and:
#'   - n columnt with the number of mutations in the VAF interval
#'   - x and y columns describing SFS
#'   - y_scaled with y values scaled to the range 0-1
#'
#' @param object SNVs tibble object
#' @param digits resolution of the cumulative tails calculation
#' @param ... other arguments
#' @examples
#' data("tcga_brca_test")
#' tcga_brca_test |>
#'   calc_SFS()
#'
#' tcga_brca_test |>
#'   plot_SFS()
#' @name sfs


#' @rdname sfs
#' @export
calc_SFS <- function(object, ...) {
  UseMethod("calc_SFS")
}


#' @rdname sfs
#' @export
plot_SFS <- function(object, ...) {
  UseMethod("plot_SFS")
}


#' Cumulative tails
#'
#' Creates  cevodata$models$cumulative_tails tbl with the groupping variables and:
#'   - n column with the number of mutations in the VAF interval
#'   - x and y columns describing the cumulative tails
#'   - y_scaled with y values scaled to the range 0-1
#'
#' @param object SNVs tibble object
#' @param digits resolution of the cumulative tails calculation
#' @param ... other arguments
#' @examples
#' data("tcga_brca_test")
#' tcga_brca_test |>
#'   calc_cumulative_tails()
#'
#' tcga_brca_test |>
#'   plot_cumulative_tails()
#' @name cumulative_tails


#' @rdname cumulative_tails
#' @export
calc_cumulative_tails <- function(object, ...) {
  UseMethod("calc_cumulative_tails")
}


#' @rdname cumulative_tails
#' @export
plot_cumulative_tails <- function(object, ...) {
  UseMethod("plot_cumulative_tails")
}
