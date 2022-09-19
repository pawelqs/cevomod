
#' @export
SNVs <- function(object, ...) {
  UseMethod("SNVs")
}


#' @export
add_SNV_data <- function(object, ...) {
  UseMethod("add_SNV_data")
}


#' @export
default_SNVs <- function(object, ...) {
  UseMethod("default_SNVs")
}


#' @export
`default_SNVs<-` <- function(object, ..., value) {
  UseMethod("default_SNVs<-")
}


#' @export
CNVs <- function(object, ...) {
  UseMethod("CNVs")
}


#' @export
add_CNV_data <- function(object, ...) {
  UseMethod("add_CNV_data")
}


#' @export
default_CNVs <- function(object) {
  UseMethod("default_CNVs")
}


#' @export
`default_CNVs<-` <- function(object, ..., value) {
  UseMethod("default_CNVs<-")
}
