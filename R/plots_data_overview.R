
#' Plot sequencing depth of cevodata
#' @param object cevodata or cevo_snvs
#' @param snvs_name name of the SNVs table
#' @param geom change geom
#' @param ... params passed to geom
#' @export
plot_sequencing_depth <- function(object, ...) {
  UseMethod("plot_sequencing_depth")
}


#' @describeIn plot_sequencing_depth Plot sequencing data of mutations in cevodata
#' @export
plot_sequencing_depth.cevodata <- function(object,
                                           snvs_name = default_SNVs(object),
                                           geom = geom_boxplot,
                                           ...) {
  snvs <- SNVs(object, snvs_name) |>
    left_join(object$metadata, by = "sample_id")
  plot_sequencing_depth(snvs, geom, ...)
}


#' @describeIn plot_sequencing_depth Plot sequencing data of mutations in cevodata
#' @export
plot_sequencing_depth.cevo_snvs <- function(object, geom = geom_boxplot, ...) {
  samples_order <- object |>
    group_by(.data$sample_id) |>
    summarise(median_dp = stats::median(.data$DP)) |>
    arrange(desc(.data$median_dp)) |>
    pull("sample_id")
  object <- object |>
    mutate(sample_id = parse_factor(.data$sample_id, levels = samples_order))
  ggplot(object, aes(.data$sample_id, .data$DP, fill = .data$sample_id)) +
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
plot_private_shared_mutations.cevodata <- function(object, geom = geom_bar, ...) {
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

  ggplot(plot_data, aes(.data$patient_id, fill = .data$n)) +
    geom(position = "fill", ...) +
    theme_minimal() +
    labs(
      y = "Fraction of mutations",
      x = "Patient",
      fill = "Samples count"
    )
}
