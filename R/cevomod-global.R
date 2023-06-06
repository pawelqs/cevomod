
cevomod_global <- new.env(parent = emptyenv())
cevomod_global$verbosity_level <- 1


#' Get verbosity level
#' @export
get_cevomod_verbosity <- function() {
  cevomod_global$verbosity_level
}


#' Change verbosity level
#' @export
set_cevomod_verbosity <- function(verbosity_level = 1) {
  old <- cevomod_global$verbosity_level
  cevomod_global$verbosity_level <- verbosity_level
  invisible(old)
}

