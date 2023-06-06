
#' Not in operator
#' @param x left-hand side argument
#' @param y right-hand side argument
#' @export
'%not in%' <- function(x,y)!('%in%'(x,y))


join_aes <- function(aes_default, aes_2) {
  aes <- c(as.list(aes_default[names(aes_default) %not in% names(aes_2)]), aes_2)
  class(aes) <- 'uneval'
  aes
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


msg <- function(...,
                collapse = "",
                col = "steelblue3",
                new_line = TRUE,
                verbose = get_cevomod_verbosity()) {
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


require_columns <- function(tbl, ...) {
  cols <- list(...) |>
    unlist()
  tbl_name <- deparse(substitute(tbl))
  missing <- cols %not in% colnames(tbl)
  if (sum(missing) > 0) {
    stop("The following columns are missing in the ", tbl_name, ": ", str_c(cols[missing], collapse = ", "))
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


#' Fill na values in the object
#' @param object object
#' @param val value to fill the NAs
#' @export
fill_na <- function(object, val) {
  object[is.na(object)] <- val
  object
}


segment <- function(vec) {
  x <- vec != lag(vec)
  x[1] <- 0
  cumsum(x)
}


chromosomes_to_int <- function(chrom) {
  case_match(
    chrom,
    "chr1" ~ 1,
    "chr2" ~ 2,
    "chr3" ~ 3,
    "chr4" ~ 4,
    "chr5" ~ 5,
    "chr6" ~ 6,
    "chr7" ~ 7,
    "chr8" ~ 8,
    "chr9" ~ 9,
    "chr10" ~ 10,
    "chr11" ~ 11,
    "chr12" ~ 12,
    "chr13" ~ 13,
    "chr14" ~ 14,
    "chr15" ~ 15,
    "chr16" ~ 16,
    "chr17" ~ 17,
    "chr18" ~ 18,
    "chr19" ~ 19,
    "chr20" ~ 20,
    "chr21" ~ 21,
    "chr22" ~ 22,
    "chrX" ~ 23,
    "chrY" ~ 24,
    "chrMT" ~ 25
  )
}


#' Quick save to ~/.cevomod directory
#' @param object object to save
#' @export
quick_save <- function(object) {
  dir.create("~/.cevomod", showWarnings = FALSE)
  write_rds(object, "~/.cevomod/object.Rds")
}


#' @describeIn quick_save Quick load of ~/.cevomod/object.Rds
#' @export
quick_load <- function() {
  read_rds("~/.cevomod/object.Rds")
}
