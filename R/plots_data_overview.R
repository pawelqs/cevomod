
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
  snvs <- SNVs(object) |>
    left_join(object$metadata, by = "sample_id")
  samples_order <- snvs |>
    group_by(.data$sample_id) |>
    summarise(median_dp = stats::median(.data$DP)) |>
    arrange(desc(.data$median_dp)) |>
    pull("sample_id")
  snvs <- snvs |>
    mutate(sample_id = parse_factor(.data$sample_id, levels = samples_order))
  default_mapping <- aes(.data$sample_id, .data$DP, fill = .data$sample_id)
  ggplot(snvs, join_aes(default_mapping, mapping)) +
    geom(...) +
    scale_y_log10() +
    theme_minimal() +
    labs(
      y = "Sequencing coverage of mutations",
      x = "Sample"
    ) +
    hide_legend()
}


#' Plot proportions of private and shared mutations
#' @param object cevodata
#' @param ... other params
#' @export
plot_private_shared_mutations <- function(object, ...) {
  UseMethod("plot_private_shared_mutations")
}


#' @describeIn plot_private_shared_mutations Plot private and shared mutation fractions
#' @inheritParams plot_sequencing_depth
#' @export
plot_private_shared_mutations.cevodata <- function(object,
                                                   mapping = NULL,
                                                   geom = geom_bar,
                                                   ...) {
  if (is.null(object$metadata$patient_id)) {
    stop("patient_id column must be present to plot this chart!")
  }
  snvs <- SNVs(object) |>
    filter(!is.na(.data$VAF), .data$VAF > 0.000001)
  plot_data <- snvs |>
    left_join(object$metadata, by = "sample_id") |>
    unite("mut_id", .data$chrom, .data$pos, .data$ref, .data$alt, sep = ":") |>
    group_by(.data$patient_id, .data$mut_id) |>
    count() |>
    mutate(n = as.character(.data$n))

  default_mapping <- aes(.data$patient_id, fill = .data$n)
  ggplot(plot_data, join_aes(default_mapping, mapping)) +
    geom(position = "fill", ...) +
    theme_minimal() +
    labs(
      y = "Fraction of mutations",
      x = "Patient",
      fill = "Samples count"
    )
}