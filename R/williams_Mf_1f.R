#' Calculate Williams M(f) ~ 1/f
#'
#' @param dt tibble with VAF column. Might be groupped and contain the data for multiple samples
#' @param digits resolution of the cumulative tails calculation
#' @return tbl with the groupping variables and:
#'   - VAF
#'   - n - number of mutations in the VAF interval
#'   - `M(f)` and `1/f` columns to plot William's statistics
#' @export
#'
#' @examples
#' data("snvs_test")
#' snvs_test |>
#'   dplyr::group_by(sample_id) |>
#'   calc_Mf_1f()
calc_Mf_1f <- function(dt, digits = 2) {
  group_variables <- group_vars(dt)
  res <- dt %>%
    mutate(VAF = round(.data$VAF, digits = digits)) %>%
    group_by(.data$VAF, .add = TRUE) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    arrange(!!!syms(group_variables), .data$VAF) %>%
    complete_missing_VAF_levels(fill = list(n = 0)) %>%
    mutate(
      `M(f)` = cumsum(n),
      `1/f` = round(1/.data$VAF, digits = 4)
    )
  class(res) <- c("cevo_Mf_1f_tbl", class(res))
  res
}


#' Plot M(f) ~ 1/f
#'
#' @param x tibble with calc_Mf_1f() results
#' @param from min VAF to plot
#' @param to max VAF to plot
#' @param scale scale data?
#' @param geom ggplot geom to use, eg. geom_line()
#' @param mapping mapping
#' @param ... futher passed to geom_()
#'
#' @return ggplot obj
#' @export
#'
#' @examples
#' data("snvs_test")
#' snvs_test |>
#'   dplyr::group_by(sample_id) |>
#'   calc_Mf_1f() |>
#'   plot()
#'
#' snvs_test |>
#'   dplyr::group_by(sample_id) |>
#'   plot_Mf_1f()
plot.cevo_Mf_1f_tbl <- function(x, from = 0.1, to = 0.25, scale = TRUE,
                                geom = "point", mapping = NULL, ...) {

  x <- filter(x, .data$VAF >= from, .data$VAF <= to)
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
    scale_x_reverse(breaks = break_vals, labels = break_labs) +
    theme_minimal() +
    labs(title = "M(f) ~ 1/f")
}


#' @describeIn plot.cevo_Mf_1f_tbl Plot M(f) ~ 1/f
#' @inherit plot.cevo_Mf_1f_tbl
#' @inherit calc_Mf_1f
#' @export
plot_Mf_1f <- function(dt, digits = 2, from = 0.1, to = 0.25, scale = TRUE, geom = "point", ...) {
  Mf_1f <- calc_Mf_1f(dt, digits)
  plot(Mf_1f, geom = geom, ...)
}
