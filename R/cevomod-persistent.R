
settings_dir <- tools::R_user_dir("cevomod", which = "config")
settings_file <- file.path(settings_dir, "settings.rds")

default_settings <- list(
  containers_dir = NULL
)


get_settings <- function() {
  if (file.exists(settings_file)) {
    readRDS(settings_file)
  } else {
    default_settings
  }
}


#' Get/Set the containers directory
#' @param dir Path for containers
#' @export
set_containers_dir <- function(dir) {
  settings <- get_settings()
  settings$containers_dir <- dir

  if(!dir.exists(settings_dir)) {
    dir.create(settings_dir)
  }
  write_rds(settings, settings_file)
}


#' @rdname set_containers_dir
#' @export
get_containers_dir <- function() {
  get_settings()$containers_dir
}
