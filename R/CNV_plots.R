
#' @describeIn cnv_plots Plot CNVs
#' @export
plot_CNV_heatmap.cevodata <- function(object, meta_field, ...) {
  granges <- CNVs(object) |>
    rename(seqnames = .data$chrom) |>
    group_by(.data$sample_id) |>
    nest() |>
    deframe() |>
    map(GenomicRanges::GRanges)

  heatmap_granges(granges, meta_field, ...)
}


# rowdata <- tibble(
#   name = names(seg_gr),
#   name2 = str_replace(name, "AMLRO_", "AMLRO"),
#   sample_short = str_replace(name2, "AMLRO_", "AMLRO")
# ) %>%
#   separate(name2, into = c("patient", "sample", "cnv"), remove = FALSE) %>%
#   unite(left, patient, sample, remove = FALSE) %>%
#   mutate(
#     patient_ = str_replace(patient, "AMLRO", "AMLRO_"),
#     sex = if_else(patient_ %in% male_patients, "male", "female")
#   )
#
# row_ha = rowAnnotation(
#   sex = rowdata$sex,
#   col = list(sex = c("male" = "#d95f02", "female" = "#7570b3"))
# )


#' @describeIn cnv_plots Plot Granges list
#' @export
heatmap_granges <- function(granges, meta_field,
                            row_groups = NULL,
                            keep_sites_present_in = floor(0.8 * length(granges)),
                            color_breaks = c(0, 2, 6),
                            colors = c("dodgerblue3", "white", "firebrick3"),
                            window_width = 1000000, upper_limit = 6,
                            cluster_rows = FALSE,
                            show_row_names = TRUE,
                            show_column_names = FALSE,
                            use_raster = TRUE,
                            cluster_columns = FALSE,
                            border = TRUE,
                            legend_params = NULL,
                            verbose = TRUE, ...) {
  score <- NULL

  ########## Get common ranges
  ranges_cov <- granges %>%
    GenomicRanges::GRangesList() %>%
    unlist() %>%
    plyranges::compute_coverage()

  range_fractions <- ranges_cov %>%
    as_tibble() %>%
    group_by(.data$score) %>%
    summarise(width = sum(.data$width)) %>%
    arrange(desc(.data$score)) %>%
    mutate(
      total_len = sum(.data$width),
      frct = .data$width / .data$total_len,
      cum = cumsum(.data$frct)
    ) %>%
    select(.data$score, .data$cum) %>%
    deframe()

  if (verbose) {
    print("Proportion of ranges present in N GRanges:")
    print(range_fractions)
    print(sprintf("Kept sites present in %d/%d GRanges", keep_sites_present_in, length(granges)))
    print(sprintf("Kept %f of all sites", range_fractions[as.character(keep_sites_present_in)]))
  }

  common_ranges <- ranges_cov %>%
    filter(score > keep_sites_present_in) %>%
    plyranges::reduce_ranges()

  ########### List of GRanges to matrix
  windows <- common_ranges %>%
    GenomicRanges::tile(width = window_width) %>%
    unlist()

  dt <- map(granges, ~ plyranges::join_overlap_left(windows, .x))
  mat <- dt %>%
    map(as_tibble) %>%
    # Summarise values for window ranges covering more than one interval
    map(~ select(.x, .data$seqnames, .data$start, .data$end, !!sym(meta_field))) %>%
    map(~ group_by(.x, .data$seqnames, .data$start, .data$end)) %>%
    map(~ summarise(.x, val = mean(!!sym(meta_field)))) %>%
    # Create matrix
    map("val") %>%
    as.data.frame() %>%
    as.matrix() %>%
    t()
  mat[mat > upper_limit] <- upper_limit

  ######### Prepare heatmap chromosome annotation
  chr_means <- as_tibble(windows) %>%
    mutate(i = row_number()) %>%
    group_by(.data$seqnames) %>%
    summarise(m = mean(.data$i)) %>%
    deframe()
  chr_labels <- vector(mode = "character", length = length(windows))
  chr_labels[chr_means] <- names(chr_means)

  seqname_replacements <- names(chr_means) %>%
    as_tibble_col(column_name = "seqname") %>%
    mutate(code = row_number() %% 2)
  chr_bin <- tibble(seqname = as.character(GenomicRanges::seqnames(windows))) %>%
    left_join(seqname_replacements, by = "seqname") %>%
    pull(.data$code)


  chr_bar <- HeatmapAnnotation(
    chr_text = ComplexHeatmap::anno_text(chr_labels, gp = gpar(fontsize = 8)),
    chr = chr_bin,
    show_legend = FALSE,
    which = "column",
    col = list(chr = c("0" = "grey88", "1" = "black"))
  )

  ########### Plot Heatmap
  Heatmap(mat,
    col = circlize::colorRamp2(breaks = color_breaks, colors),
    row_split = if (is.null(row_groups)) NULL else factor(row_groups, levels = unique(row_groups)),
    cluster_row_slices = FALSE,
    top_annotation = chr_bar,
    heatmap_legend_param = c(list(title = meta_field), legend_params),
    cluster_rows = cluster_rows,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    use_raster = use_raster,
    cluster_columns = cluster_columns,
    border = border,
    ...
  )
}
