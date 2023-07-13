
cevomod_global <- new.env(parent = emptyenv())
cevomod_global$verbosity_level <- 1


#' Get the verbosity level
#' @export
get_cevomod_verbosity <- function() {
  cevomod_global$verbosity_level
}


#' Change the verbosity level
#' @param verbosity_level Verbosity level to use:
#'   0 - silent
#'   1 - normal
#'   2 - detailed (in some cases)
#' @export
set_cevomod_verbosity <- function(verbosity_level = 1) {
  old <- cevomod_global$verbosity_level
  cevomod_global$verbosity_level <- verbosity_level
  invisible(old)
}

