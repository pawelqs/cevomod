
#' Plot sequencing depth of cevodata
#' @param object cevodata
#' @param ... other params
#' @export
plot_sequencing_depth <- function(object, ...) {
  UseMethod("plot_sequencing_depth")
}


#' @describeIn plot_sequencing_depth Plot sequencing data of mutations in cevodata
#' @param mapping customize mapping
#' @param geom change geom
#' @param ... params passed to geom
#' @export
plot_sequencing_depth.cevodata <- function(object,
                                           mapping = NULL,
                                           geom = geom_boxplot,
                                           ...) {
  snvs <- SNVs(object)
  default_mapping <- aes(.data$sample_id, .data$DP, fill = .data$sample_id)
  ggplot(snvs, join_aes(default_mapping, mapping)) +
    geom(...) +
    scale_y_log10() +
    theme_minimal() +
    labs(
      y = "Sequencing depth of mutations",
      x = "Sample"
    ) +
    hide_legend()
}
