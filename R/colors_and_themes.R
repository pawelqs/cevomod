
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
    select("pal_name", "colors") %>%
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


# ----------------------------- Theme presettings ----------------------------

#' Hide legend
#' @export
hide_legend <- function() {
  theme(legend.position = "none")
}


#' Rotate x axix labels
#' @param angle angle
#' @export
rotate_x_labels <- function(angle = 90) {
  theme(axis.text.x = element_text(angle = angle))
}


#' Hide facet labels
#' @export
hide_facet_labels <- function() {
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
}

# ---------------------------- color scales ------------------------------------

colors <- list(
  starfleet = c("#5B1414", "#AD722C", "#1A6384")
)


# ggplot2::discrete_scale("color", "ellie", )


#' Different color pallettes
#' @param palette palette name
#' @param ... other arguments passed to scale_*_manual()
#' @export
scale_fill_cevomod <- function(palette = "starfleet", ...) {
  scale_fill_manual(values = colors[[palette]], ...)
}


#' @rdname scale_fill_cevomod
#' @export
scale_color_cevomod <- function(palette = "starfleet", ...) {
  scale_color_manual(values = colors[[palette]], ...)
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
