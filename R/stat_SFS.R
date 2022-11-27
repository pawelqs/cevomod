
#' Site Frequency Spectra
#'
#' Creates  cevodata$models$SFS with the groupping variables and:
#'   - n columnt with the number of mutations in the VAF interval
#'   - x and y columns describing SFS
#'   - y_scaled with y values scaled to the range 0-1
#'
#' @param object SNVs tibble object
#' @param digits resolution of the cumulative tails calculation
#' @param ... other arguments
#' @examples
#' data("tcga_brca_test")
#' tcga_brca_test |>
#'   calc_SFS()
#'
#' tcga_brca_test |>
#'   plot_SFS() +
#'   layer_mutations(drivers = "BRCA")
#'
#' SNVs(tcga_brca_test) |>
#'   dplyr::group_by(sample_id) |>
#'   calc_SFS() |>
#'   plot()
#' @name sfs


#' @rdname sfs
#' @export
calc_SFS <- function(object, ...) {
  UseMethod("calc_SFS")
}


#' @describeIn sfs Calculate SFS
#' @export
calc_SFS.cevodata <- function(object, digits = 2, ...) {
  object$models[["SFS"]] <- SNVs(object) |>
    group_by(.data$sample_id) |>
    calc_SFS(digits = 2) |>
    ungroup()
  object
}


#' @describeIn sfs Calculate SFS
#' @export
calc_SFS.tbl_df <- function(object, digits = 2, ...) {
  group_variables <- group_vars(object)
  res <- object %>%
    mutate(VAF = round(.data$VAF, digits = digits)) %>%
    group_by(.data$VAF, .add = TRUE) %>%
    summarise(y = n(), .groups = "drop_last") %>%
    arrange(.data$VAF, .by_group = TRUE) %>%
    complete_missing_VAF_levels(fill = list(y = 0)) %>%
    mutate(y_scaled = round(.data$y / sum(.data$y), digits = 4))
  class(res) <- c("cevo_SFS_tbl", class(res))
  res
}


#' @rdname sfs
#' @export
plot_SFS <- function(object, ...) {
  UseMethod("plot_SFS")
}


#' Plot SFS
#'
#' @param x tibble with calc_SFS() results
#' @param y_scaled logical
#' @param ... futher passed to geom_()
#' @return ggplot obj
#' @export
plot.cevo_SFS_tbl <- function(x, y_scaled = FALSE, ...) {
  group_variables <- group_vars(x)
  y <- if (y_scaled) "y_scaled" else "y"
  n_colors <- n_distinct(x[group_variables[[1]]])

  ggplot(x) +
    aes(.data$VAF, !!sym(y), color = !!sym(group_variables[[1]]), ...) +
    geom_line() +
    theme_ellie(n_colors) +
    labs(
      title = "SFS",
      y = "count"
    )
}


#' @describeIn sfs Plot SFS
#' @export
plot_SFS.cevodata <- function(object, ..., geom = "line", alpha = if (geom == "bar") 0.8 else 1) {
  dt <- SNVs(object) |>
    filter(.data$VAF > 0.00001) |>
    left_join(object$metadata, by = "sample_id")
  ggplot(dt, aes(.data$VAF, color = .data$sample_id, fill = .data$sample_id)) +
    stat_SFS(..., geom = geom, alpha = alpha) +
    theme_ellie(n = n_distinct(dt$sample_id)) +
    labs(
      title = "SFS",
      y = "count"
    )
}


#' SFS stat
#'
#' Modification of stat_bin with custom parameters
#' @param mapping aes()
#' @param data data
#' @param binwidth binwidth
#' @param geom geom
#' @param position position
#' @param ... args passed to stat_bin
#'
#' @export
stat_SFS <- function(mapping = NULL, data = NULL,
                     binwidth = 0.01, geom = "line", position = "identity", ...) {
  stat_bin(mapping = mapping, data = data, binwidth = binwidth, geom = geom, position = position, ...)
}
