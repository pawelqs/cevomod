
settings_dir <- tools::R_user_dir("cevomod", which = "config")
settings_file <- file.path(settings_dir, "settings.yml")


get_settings <- function() {
  if (file.exists(settings_file)) {
    yaml::read_yaml(settings_file)
  } else {
    list(
      containers_dir = NULL
    )
  }
}


#' Set the containers directory
#' @param dir Path for containers
#' @export
set_containers_dir <- function(dir) {
  settings <- get_settings()
  settings$containers_dir <- dir

  if(!dir.exists(settings_dir)) {
    dir.create(settings_dir)
  }
  yaml::write_yaml(settings, settings_file)
}


#' Get the verbosity level
#' @export
get_containers_dir <- function() {
  get_settings()$containers_dir
}
