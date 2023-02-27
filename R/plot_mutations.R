
#' Show mutations in particular genes
#'
#' @param object cevodata
#' @param genes list of genes for which mutations should be shown
#' @param drivers cancer type, for which mutations in driver genes should be shown. Needs to
#'   be taken from `driver_genes`
#' @param mark_genes list of genes to mark
#' @param y which to show on y axis: "genes" or "samples"
#' @param shape "impact" or "variant_classification"?
#' @param filter_fun Function for filtering mutations:
#'   variant_classification_filter() or impact_filter()
#' @return ggplot obj
#' @name mutation_plots


#' @rdname mutation_plots
#' @export
plot_mutations <- function(object, ...) {
  UseMethod("plot_mutations")
}


#' @rdname mutation_plots
#' @export
plot_mutations.cevodata <- function(object, genes = NULL, drivers = NULL, mark_genes = NULL,
                           y = "genes", shape = "impact",
                           filter_fun = guess_filter_fun(shape), ...) {
  genes_data <- SNVs(object) %>%
    filter_SNVs(genes, drivers)

  if (!is.null(filter_fun)) {
    genes_data <- filter_fun(genes_data)
  }

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


#' @rdname mutation_plots
#' @export
plot_mutations.tbl_df <- function(object, genes = NULL, drivers = NULL, mark_genes = NULL,
                                  y = "genes", shape = "impact",
                                  filter_fun = guess_filter_fun(shape), ...) {
  genes_data <- object %>%
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


#' @describeIn mutation_plots Adds mutations to SFS plots
#' @param color color
#' @param size size
#' @param show_labels `lgl` use ggrepel to label the mutations?
#' @param ... other arguments passed to geom_point()
#' @export
layer_mutations <- function(object,
                            genes = NULL, drivers = NULL, show_labels = TRUE,
                            color = "black", size = 3, shape = "impact",
                            filter_fun = guess_filter_fun(shape), ...) {
  dt <- SNVs(object)
  list(
    scale_shape_manual(values = c(2, 4, 3, 5)),
    geom_point(
      aes(x = .data$VAF, shape = .data[[shape]]),
      data = dt %>%
        filter_SNVs(genes, drivers) %>%
        filter_fun(),
      y = 0,
      size = size,
      color = color,
      inherit.aes = FALSE,
      ...
    ),
    if (show_labels) {
      rlang::check_installed("ggrepel", reason = "to label driver mutations on the plot")
      ggrepel::geom_label_repel(
        aes(x = .data$VAF, shape = .data[[shape]], label = .data$gene_symbol),
        data = dt %>%
          filter_SNVs(genes = NULL, drivers) %>%
          filter_fun(),
        y = 0,
        size = size,
        inherit.aes = FALSE,
        max.overlaps = 15
      )
    }
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
