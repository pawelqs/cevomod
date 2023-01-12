'%not in%' <- function(x,y)!('%in%'(x,y))


join_aes <- function(aes_default, aes_2) {
  aes <- c(as.list(aes_default[names(aes_default) %not in% names(aes_2)]), aes_2)
  class(aes) <- 'uneval'
  aes
}


complete_missing_VAF_levels <- function(dt, fill, digits = 2) {
  VAF <- NULL
  group_variables <- group_vars(dt)
  dt %>%
    mutate(VAF = round(VAF, digits = digits)) %>%
    mutate(
      VAF = VAF %>%
        as.character() %>%
        parse_factor(levels = as.character(seq(0, 1, by = 1/(10^digits))))
    ) %>%
    complete(VAF, fill = fill) %>%
    mutate(
      VAF = VAF %>%
        as.character() %>%
        parse_double()
    ) %>%
    group_by(!!!syms(group_variables))
}


drop_na_columns <- function(.data) {
  .data |>
    keep(~all(!is.na(.x)))
}


get_VAF_range <- function(snvs, pct_left = 0.05, pct_right = 0.95) {
  bounds <- snvs |>
    filter(.data$VAF > 0.00001, !is.na(.data$VAF)) |>
    group_by(.data$sample_id) |>
    summarise(
      lower_bound = stats::quantile(.data$VAF, pct_left),
      higher_bound = stats::quantile(.data$VAF, pct_right)
    )
  bounds
}


msg <- function(..., collapse = "", col = "steelblue3", new_line = TRUE, verbose = TRUE) {
  msg <- str_c(list(...), collapse = collapse)
  if (verbose && new_line) {
    cli::cat_line(msg, col = col)
  } else if (verbose) {
   cat(crayon::blue(msg))
  }
}


require_packages <- function(...) {
  pkgs <- list(...)
  missing <- !map_lgl(pkgs, requireNamespace, quietly = TRUE)
  if (any(missing)) {
    stop(
      paste0("Package '", pkgs[missing], "' must be installed to use this function.\n"),
      call. = FALSE
    )
  }
}


#' Run cevobrowser app
#' @export
run_browser <- function() {
  require_packages("shiny", "shinydashboard", "shinyWidgets")

  app_dir <- system.file("cevobrowser", package = "cevomod")
  if (app_dir == "") {
    stop("Could not find app directory. Try re-installing `cevomod`.", call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal")
}


#' @export
print.cevo_snvs <- function(x, ...) {
  msg("<cevo_snvs> tibble")
  NextMethod()
}


#' Shuffle order of elements in object
#' @param object object to shuffle
#' @param ... other arguments
#' @export
shuffle <- function(object, ...) {
  UseMethod("shuffle")
}


#' @describeIn shuffle Shuffle order of rows in tibble
#' @export
#' @examples
#' tibble::tibble(i = 1:10) |>
#'   shuffle()
shuffle.tbl_df <- function(object, ...)  {
  object[sample(1:nrow(object)), ]
}
