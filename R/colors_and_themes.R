
# ---------------------------- color scales ------------------------------------


palettes <- structure(
  list(
    starfleet = c("#5B1414", "#AD722C", "#1A6384")
  ),
  class = c("cevo_palettes", "list")
)


#' Show implemented palettes
#' @param x Palettes
#' @param ... other args
#' @export
print.cevo_palettes <- function(x, ...) {
  x |>
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



#' Different color pallettes
#' @param palette palette name
#' @param ... other arguments passed to scale_*_manual()
#' @export
scale_fill_cevomod <- function(palette = "starfleet", ...) {
  scale_fill_manual(values = palettes[[palette]], ...)
}


#' @rdname scale_fill_cevomod
#' @export
scale_color_cevomod <- function(palette = "starfleet", ...) {
  scale_color_manual(values = palettes[[palette]], ...)
}


# ----------------------------- Theme presettings ----------------------------

#' Hide legend
#' @export
hide_legend <- function() {
  theme(legend.position = "none")
}


#' Rotate x axix labels
#' @param angle angle
#' @param vjust vjust
#' @export
rotate_x_labels <- function(angle = 90, vjust = 0.5) {
  theme(axis.text.x = element_text(angle = 90, vjust = vjust))
}


#' Hide facet labels
#' @export
hide_facet_labels <- function() {
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
}
