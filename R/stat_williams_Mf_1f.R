
#' William's M(f) ~ 1/f statistics
#'
#' Mf_1f columns description:
#' - sample_id
#' - f
#' - n - number of mutations in the VAF interval
#' - `M(f)` and `1/f` columns to plot William's statistics
#'
#' @param object SNVs tibble object
#' @param which_snvs Which SNVs to use?
#' @param column VAF or CCF/2
#' @param bins Resolution of the cumulative tails calculation
#' @param verbose Verbose?
#' @param ... Other arguments
#' @examples
#' data("test_data")
#' test_data |>
#'   calc_Mf_1f()
#'
#' test_data |>
#'   plot_Mf_1f()
#' @name Mf_1f
NULL


#' @describeIn Mf_1f Calculate Williams M(f) ~ 1/f and saves it to
#'   cevodata$models$Mf_1f
#' @export
calc_Mf_1f <- function(object, ...) {
  UseMethod("calc_Mf_1f")
}


#' @describeIn Mf_1f Method for <cevodata> object
#' @export
calc_Mf_1f.cevodata <- function(object,
                                which_snvs = default_SNVs(object),
                                column = get_frequency_measure_name(object, which_snvs),
                                bins = 100,
                                verbose = get_cevomod_verbosity(),
                                ...) {
  Mf_1f <- SNVs(object, which_snvs) |>
    calc_Mf_1f(column = column, bins = bins, verbose = verbose)
  object$models[["Mf_1f"]] <- Mf_1f
  object
}


#' @describeIn Mf_1f Method for <cevo_snvs> object
#' @export
calc_Mf_1f.cevo_snvs <- function(object,
                                 column = get_frequency_measure_name(object),
                                 bins = 100,
                                 verbose = get_cevomod_verbosity(),
                                 ...) {
  msg("Calculating Williams's M(f) ~ 1/f statistics, using ", column, " column", verbose = verbose)

  snvs <- cut_f_intervals(object, bins = bins, column = column)
  intervals <- attributes(snvs)$intervals
  res <- snvs |>
    group_by(.data$sample_id, .data$f_interval) %>%
    summarise(n = n(), .groups = "drop_last") |>
    complete_missing_f_intervals(intervals) |>
    replace_na(list(n = 0)) |>
    mutate(f = get_interval_centers(.data$f_interval), .after = "f_interval") |>
    arrange(desc(.data$f), .by_group = TRUE) %>%
    mutate(
      `M(f)` = cumsum(n),
      `1/f` = round(1/.data$f, digits = 4)
    ) |>
    ungroup()
  class(res) <- c("cevo_Mf_1f_tbl", class(res))
  res
}


#' @rdname Mf_1f
#' @export
plot_Mf_1f <- function(object, ...) {
  UseMethod("plot_Mf_1f")
}


#' @describeIn Mf_1f Plot M(f) ~ 1/f
#' @inherit plot.cevo_Mf_1f_tbl
#' @export
plot_Mf_1f.cevodata <- function(object,
                                bins = NULL, from = 0.1, to = 0.25,
                                scale = TRUE, geom = "point", ...) {
  get_Mf_1f(object) |>
    left_join(object$metadata, by = "sample_id") |>
    plot(geom = geom, from = from, to = to, scale = scale, ...)
}


#' Plot M(f) ~ 1/f
#'
#' @param x tibble with calc_Mf_1f() results
#' @param from min f to plot
#' @param to max f to plot
#' @param scale scale data?
#' @param geom ggplot geom to use, eg. geom_line()
#' @param mapping mapping
#' @param ... futher passed to geom_()
#'
#' @return ggplot obj
#' @export
#'
#' @examples
#' data("test_data")
#'
#' test_data |>
#'   calc_Mf_1f() |>
#'   plot_Mf_1f()
plot.cevo_Mf_1f_tbl <- function(x, from = 0.1, to = 0.25, scale = TRUE,
                                geom = "point", mapping = NULL, ...) {

  x <- x |>
    filter(.data$f >= from, .data$f <= to) |>
    group_by(.data$sample_id)
  if (scale) {
    x <- x %>%
      mutate(
        `M(f)` = .data$`M(f)` - min(.data$`M(f)`),
        `M(f)` = .data$`M(f)` / max(.data$`M(f)`),
      )
  }

  breaks <- seq(from, to, length.out = 4) %>%
    round(digits = 2) %>%
    rev()
  break_vals <- 1 / breaks
  break_labs  <- str_c("1/", breaks)
  default_aes <- aes(.data$`1/f`, .data$`M(f)`, color = .data$sample_id)
  mapping <- join_aes(default_aes, mapping)

  geom <- list(
    mapping,
    if (geom == "line") geom_line(...),
    if (geom == "point") geom_point(size = 1, ...)
  )

  ggplot(x) +
    geom +
    scale_x_continuous(breaks = break_vals, labels = break_labs) +
    theme_minimal() +
    labs(title = "M(f) ~ 1/f")
}


#' @describeIn Mf_1f Get Mf_1f
#' @param model_name name of slot with Mf_1f statistics
#' @export
get_Mf_1f <- function(object, model_name = "Mf_1f", verbose = TRUE, ...) {
  Mf_1f <- object$models[[model_name]]
  if (is.null(Mf_1f)) {
    msg(
      "Mf_1f's not calculated yet. Calculating with default bins",
      verbose = verbose
    )
    object <- calc_Mf_1f(object)
    Mf_1f <- object$models[[model_name]]
  }
  Mf_1f
}
