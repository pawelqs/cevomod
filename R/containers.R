

#' Build the CliP Apptainer container
#'
#' CliP.sif is saved to
#'   - out_dir, if provided
#'   - if out_dir is NULL but the containers_dir was set using the [set_containers_dir()],
#'     image will be saved to the set containers_dir.
#'   - if out_dir is NULL and  the containers_dir was not set, image will
#'     be saved to the current working dir.
#' @param out_dir Path
#' @export
build_clip_container <- function(out_dir = NULL, force = FALSE) {
  if (!is_apptainer_installed()) {
    stop("Apptainer needs to be installed to build the containers")
  }

  if (is.null(out_dir) & !is.null(get_containers_dir())) {
    out_dir <- get_containers_dir()
  } else {
    out_dir <- "."
  }

  sif_file <- file.path(out_dir, "CliP.sif") |>
    str_replace("//", "/")

  if (!file.exists(sif_file) | force) {
    def_file <- system.file("CliP.def", package = "cevomod")
    command <- str_c("apptainer build", sif_file, def_file, sep = " ")
    system(command)
  } else {
    msg(sif_file, " exists. Set force = TRUE to overwrite it")
  }
}


is_apptainer_installed <- function() {
  rlang::check_installed("processx", reason = "to interact with system")
  x <- processx::run("apptainer", args = "--version")
  x$status == 0
}



