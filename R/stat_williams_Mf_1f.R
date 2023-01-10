
#' William's M(f) ~ 1/f statistics
#'
#' Creates  cevodata$models$Mf_1f tbl with the groupping variables and:
#'   - VAF
#'   - n - number of mutations in the VAF interval
#'   - `M(f)` and `1/f` columns to plot William's statistics
#'
#' @param object SNVs tibble object
#' @param bins resolution of the cumulative tails calculation
#' @param ... other arguments
#' @examples
#' data("tcga_brca_test")
#' tcga_brca_test |>
#'   calc_Mf_1f()
#'
#' tcga_brca_test |>
#'   plot_Mf_1f()
#' @name Mf_1f
NULL


#' @describeIn Mf_1f Calculate Williams M(f) ~ 1/f
#' @export
calc_Mf_1f <- function(object, ...) {
  UseMethod("calc_Mf_1f")
}


#' @describeIn Mf_1f Calculate Williams M(f) ~ 1/f
#' @export
calc_Mf_1f.cevodata <- function(object, bins = 100, ...) {
  Mf_1f <- SNVs(object) |>
    calc_Mf_1f(bins = bins)
  # class(Mf_1f) <- c("cevo_Mf_1f_tbl", class(Mf_1f))
  object$models[["Mf_1f"]] <- Mf_1f
  object
}


#' @describeIn Mf_1f Calculate Williams M(f) ~ 1/f
#' @export
calc_Mf_1f.cevo_snvs <- function(object, bins = NULL, ...) {
  snvs <- cut_VAF_intervals(object, bins = bins)
  intervals <- attributes(snvs)$intervals
  res <- snvs |>
    group_by(.data$sample_id, .data$VAF_interval) %>%
    summarise(n = n(), .groups = "drop_last") |>
    complete_missing_VAF_intervals(intervals) |>
    replace_na(list(n = 0)) |>
    mutate(VAF = get_interval_centers(.data$VAF_interval), .after = "VAF_interval") |>
    arrange(desc(.data$VAF), .by_group = TRUE) %>%
    mutate(
      `M(f)` = cumsum(n),
      `1/f` = round(1/.data$VAF, digits = 4)
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
#' data("tcga_brca_test")
#'
#' tcga_brca_test |>
#'   calc_Mf_1f() |>
#'   plot_Mf_1f()
plot.cevo_Mf_1f_tbl <- function(x, from = 0.1, to = 0.25, scale = TRUE,
                                geom = "point", mapping = NULL, ...) {

  x <- x |>
    filter(.data$VAF >= from, .data$VAF <= to) |>
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


#' @describeIn Mf_1f Plot M(f) ~ 1/f
#' @inherit plot.cevo_Mf_1f_tbl
#' @export
plot_Mf_1f <- function(object, bins = NULL, from = 0.1, to = 0.25, scale = TRUE, geom = "point", ...) {
  plot(object$models$Mf_1f, geom = geom, from = from, to = to, scale = scale, ...)
}
