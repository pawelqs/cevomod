#' Show mutations in particular genes
#'
#' @param data mutations tibble
#' @param genes list of genes for which mutations should be shown
#' @param drivers cancer type, for which mutations in driver genes should be shown. Needs to
#'   be taken from `driver_genes`
#' @param mark_genes list of genes to mark
#' @param y which to show on y axis: "genes" or "samples"
#' @param shape "impact" or "variant_classification"?
#' @param filter_fun Function for filtering mutations:
#'   variant_classification_filter() or impact_filter()
#'
#' @return ggplot obj
#' @export
plot_mutations <- function(data, genes = NULL, drivers = NULL, mark_genes = NULL,
                           y = "genes", shape = "impact",
                           filter_fun = guess_filter_fun(shape)) {
  genes_data <- data %>%
    filter_SNVs(genes, drivers) %>%
    filter_fun()

  mark_genes_data <-
    if (!is.null(mark_genes)) {
      filter(genes_data, .data$gene_symbol %in% mark_genes)
    } else {
      filter(genes_data, FALSE)
    }

  elements <- list(
    if (y == "genes") aes(y = .data$gene_symbol, color = .data$sample_id),
    if (y == "genes") {
      scale_color_manual(
        values = pnw_palette(name = "Starfish", n_distinct(genes_data$sample_id), type = "continuous")
      )
    },
    if (y == "samples") aes(y = .data$sample_id),
    if (nrow(mark_genes_data) > 0) geom_point(aes(color = .data$gene_symbol), data = mark_genes_data),
    if (nrow(mark_genes_data) > 0) scale_color_brewer(palette = "Dark2")
  )

  ggplot(genes_data, aes(x = .data$VAF, shape = .data[[shape]])) +
    geom_point() +
    elements +
    xlim(c(0, 1)) +
    theme_minimal()
}


filter_SNVs <- function(dt, genes = NULL, drivers = NULL) {
  if (!is.null(drivers)) {
    genes <- driver_genes %>%
      filter(.data$Cancer == drivers) %>%
      pull(.data$Gene) %>%
      c(genes)
  }
  filter(dt, .data$gene_symbol %in% genes)
}


#' Adds mutations to SFS plots
#'
#' @inherit plot_mutations
#' @param color color
#' @param size size
#' @param ... other arguments passed to geom_point()
#' @export
geom_mutations <- function(genes = NULL, drivers = NULL, data = NULL,
                           color = "black", size = 3, shape = "impact",
                           filter_fun = guess_filter_fun(shape), ...) {
  . <- NULL
  list(
    scale_shape_manual(values = c(2, 4, 3, 5)),
    geom_point(
      aes(x = .data$VAF, shape = .data[[shape]]),
      data = . %>%
        filter_SNVs(genes, drivers) %>%
        filter_fun(),
      y = 0,
      size = size,
      color = color,
      ...
    )
  )
}


guess_filter_fun <- function(shape) {
  if (shape == "impact") {
    impact_filter()
  } else if (shape == "variant_classification") {
    variant_classification_filter()
  } else {
    function(dt) filter(dt, TRUE)
  }
}


#' Filter mutations by variant_classification field
#' @param keep values to keep
#' @return filtering function
#' @export
variant_classification_filter <- function(keep = c("nonsilent", "null")) {
  function(dt) {
    filter(dt, .data$variant_classification %in% keep)
  }
}


#' @describeIn variant_classification_filter Filter mutations by impact field
#' @export
impact_filter <- function(keep = c("HIGH", "MODERATE")) {
  function(dt) {
    filter(dt, .data$impact %in% keep)
  }
}
