
#' Intratumor heterogeneity
#' @param object object
#' @param ... not used
#' @name ITH
NULL


#' @rdname ITH
#' @export
plot_ITH <- function(object, ...) {
  UseMethod("plot_ITH")
}


#' @rdname ITH
#' @export
plot_ITH.cevodata <- function(object, ...) {
  ITH <- estimate_ITH(object)
  patient_metadata <- get_patients_data(object$metadata)
  ITH |>
    left_join(patient_metadata, by = "patient_id") |>
    plot()
}


#' @rdname ITH
#' @param x ITH tibble to plot
#' @export
plot.cevo_ITH_tbl <- function(x, ...) {
  x |>
    arrange(desc(.data$Jaccard_index)) |>
    mutate(patient_id = parse_factor(.data$patient_id)) |>
    ggplot(aes(.data$patient_id, .data$Jaccard_index)) +
    geom_point() +
    coord_cartesian(ylim = c(0, 1)) +
    theme(axis.text.x = element_text(angle = 90))
}


#' @rdname ITH
#' @export
estimate_ITH <- function(object, ...) {
  UseMethod("estimate_ITH")
}


#' @export
estimate_ITH.cevodata <- function(object, ...) {
  dt <- SNVs(object) |>
    unite_mutation_id() |>
    filter(.data$alt_reads > 0) |>
    left_join(object$metadata, by = "sample_id") |>
    select("patient_id", "sample", "mutation_id") |>
    mutate(sample = as.character(.data$sample)) |>
    nest_by(.data$patient_id)
  ITH <- dt |>
    summarise(calc_Jaccard_index(.data$data), .groups = "drop") |>
    rename(sample1 = "group1", sample2 = "group2")
  structure(ITH, class = c("cevo_ITH_tbl", class(ITH)))
}


calc_Jaccard_index <- function(tbl) {
  groups_vec <- tbl[[1]]
  items_vec <- tbl[[2]]
  groups <- unique(tbl[[1]])
  res <- get_combinations_tbl(groups) |>
    set_names(c("group1", "group2")) |>
    mutate(Jaccard_index = NA_real_)
  for(i in 1:nrow(res)) {
    A <- items_vec[groups_vec == res$group1[[i]]]
    B <- items_vec[groups_vec == res$group2[[i]]]
    res$Jaccard_index[[i]] <- length(intersect(A, B)) / length(union(A, B))
  }
  res
}


get_combinations_tbl <- function(items) {
  utils::combn(items, m = 2) |>
    t() |>
    `colnames<-`(c("item1", "item2")) |>
    as_tibble()
}
