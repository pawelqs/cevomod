
#' @describeIn Mf_1f Calculate Williams M(f) ~ 1/f
#' @export
calc_Mf_1f.cevodata <- function(object, digits = 2, ...) {
  Mf_1f <- SNVs(object) |>
    group_by(.data$sample_id) |>
    calc_Mf_1f() |>
    ungroup()
  class(Mf_1f) <- c("cevo_Mf_1f_tbl", class(Mf_1f))
  object$models[["Mf_1f"]] <- Mf_1f
  object
}


#' @describeIn Mf_1f Calculate Williams M(f) ~ 1/f
#' @export
calc_Mf_1f.tbl_df <- function(object, digits = 2, ...) {
  group_variables <- group_vars(object)
  res <- object %>%
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
#' data("tcga_brca_test")
#'
#' SNVs(tcga_brca_test) |>
#'   dplyr::group_by(sample_id) |>
#'   calc_Mf_1f() |>
#'   plot()
#'
#' plot_Mf_1f(tcga_brca_test)
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


#' @describeIn Mf_1f Plot M(f) ~ 1/f
#' @inherit plot.cevo_Mf_1f_tbl
#' @export
plot_Mf_1f <- function(object, digits = 2, from = 0.1, to = 0.25, scale = TRUE, geom = "point", ...) {
  Mf_1f <- SNVs(object) |>
    group_by(.data$sample_id) |>
    calc_Mf_1f(digits) |>
    left_join(object$metadata, by = "sample_id") |>
    group_by(.data$patient_id, .data$sample_id, .data$sample)
  class(Mf_1f) <- c("cevo_Mf_1f_tbl", class(Mf_1f))
  plot(Mf_1f, geom = geom, from = from, to = to, scale = scale,  ...)
}
