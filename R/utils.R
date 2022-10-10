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
      lower_bound = quantile(VAF, pct_left),
      higher_bound = quantile(VAF, pct_right)
    )
  bounds
}


#' Run cevobrowser app
#' @export
run_browser <- function() {
  app_dir <- system.file("cevobrowser", package = "cevomod")
  if (app_dir == "") {
    stop("Could not find app directory. Try re-installing `cevomod`.", call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal")
}
