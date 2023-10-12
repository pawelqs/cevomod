

#' Build Apptainer container for CliP
#' @param out_dir path to save the CliP.sif file
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
  system("apptainer --version") == 0
}



