
#' Build Apptainer container for CliP
#' @param out_dir path to save the CliP.sif file
#' @export
build_CliP_container <- function(out_dir = ".") {
  sif_file <- file.path(out_dir, "CliP.sif")
  def_file <- system.file("CliP.def", package = "cevomod")
  command <- str_c("apptainer build", sif_file, def_file, sep = " ")
  system(command)
}

