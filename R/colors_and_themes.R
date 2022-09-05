
#' @importFrom paletteer scale_color_paletteer_d scale_fill_paletteer_d paletteer_d
#' @importFrom paletteer scale_color_paletteer_c scale_fill_paletteer_c paletteer_c
NULL


#' Show available palettes from PNWColorspackage
#' @return palette names
#' @export
list_pnw_palettes <- function() {
  names(PNWColors::pnw_palettes)
}


#' Show available palettes from paletteer
#' @param packages list of packages
#' @return tbl with palettes
#' @export
list_discrete_palettes <- function(packages = c("PNWColors")) {
  if (is.null(packages))
    palettes_d_names
  else
    palettes_d_names %>%
    filter(.data$package %in% packages)
}


#' @inherit list_discrete_palettes
#' @export
list_continuous_palettes <- function(packages = c("gameofthrones", "harrypotter")) {
  if (is.null(packages))
    paletteer::palettes_c_names
  else
    paletteer::palettes_c_names %>%
    filter(.data$package %in% packages)
}


#' Use PNWColors palette
#' @param palette palette name
#' @param direction 1/-1
#' @param dynamic TRUE/FALSE
#' @param ... additional arguments to pass to discrete_scale
#' @export
scale_color_pnw <- function(palette = "Sunset", direction = 1,
                           dynamic = FALSE, ...) {
  palette <- str_c("PNWColors::", palette)
  scale_color_paletteer_d(palette, direction, dynamic, ...)
}


#' @inherit scale_color_pnw
#' @export
scale_fill_pnw <- function(palette = "Sunset", direction = 1,
                            dynamic = FALSE, ...) {
  palette <- str_c("PNWColors::", palette)
  scale_fill_paletteer_d(palette, direction, dynamic, ...)
}


#' Show palettes
#'
#' @param packages list of packages
#' @param n number of colors
#' @export
print_palettes <- function(packages = c("PNWColors", "nord"), n = 10) {
  palettes <- list(
      discrete = palettes_d_names,
      continuous = palettes_c_names,
      dynamic = palettes_dynamic_names
    ) %>%
    bind_rows(.id = "dc") %>%
    filter(.data$package %in% packages) %>%
    mutate(pal_name = str_c(.data$package, .data$palette, sep = "::"))

  palettes$colors <- pmap(palettes, function(pal_name, dc, ...) {
    if (dc == "discrete")
      paletteer_d(pal_name)
    else
      paletteer_d(pal_name, n)
  })

  palettes %>%
    select(.data$pal_name, .data$colors) %>%
    deframe() %>%
    iwalk(pretty_color_print, pad_width = 25)
}


pretty_color_print <- function(x, pal_name, pad_width) {
  # Function from prismatic
  cols <- vapply(x, color_styler, FUN.VALUE = character(1), USE.NAMES = FALSE)
  cat(str_pad(pal_name, width = pad_width, side = "right"))
  cat(paste(c(cols, "\n"), collapse = " "))
}


color_styler <- function(x) {
  # Function from prismatic
  text <- crayon::make_style(prismatic::best_contrast(x), bg = FALSE)
  background <- crayon::make_style(x, bg = TRUE, colors = 256, grey = FALSE)
  crayon::combine_styles(text, background)(x)
}


#' Theme Ellie
#' @param n number of colors
#' @export
#' @importFrom PNWColors pnw_palette
theme_ellie <- function(n = 6) {
  n <- if (n < 5) 5 else n
  list(
    theme_minimal(),
    # scale_color_paletteer_d("Starfish", direction = -1),
    # scale_fill_paletteer_d("Starfish", direction = -1)
    scale_fill_manual(values = pnw_palette(name = "Starfish", n, type = "continuous")),
    scale_color_manual(values = pnw_palette(name = "Starfish", n, type = "continuous"))
  )
}


# ggplot2::discrete_scale("color", "ellie", )

#' Remove color scale
#' @export
hide_legend <- function() {
  theme(legend.position = "none")
}
