
#' @export
dplyr::filter

#' Filter/subset cevodata object
#'
#' This is a wrapper around dplyr::filter function which can be used to subset
#' cevodata object. Works like dplyr::filter, performs the filtering on metadata,
#' then filters SNVs, CNVs, clones and models keeping samples kept in metadata
#'
#' @param .data cevodata object
#' @param ... expression passed to dplyr::filter(...)
#' @inheritParams dplyr::filter
#' @return cevodata object
#' @export
filter.cevodata <- function(.data, ..., .preserve = FALSE) {
  new_object <- .data
  new_object$metadata <- filter(new_object$metadata, ...)
  ids <- new_object$metadata$sample_id
  new_object$SNVs <- map(new_object$SNVs, ~filter(.x, sample_id %in% ids))
  new_object$CNVs <- map(new_object$CNVs, ~filter(.x, sample_id %in% ids))
  new_object$clones <- map(new_object$clones, ~filter(.x, sample_id %in% ids))
  new_object$models <- map(new_object$models, ~filter(.x, sample_id %in% ids))
  new_object
}


#' Merge two cevodata objects
#' @inheritParams base::merge
#' @param name Name of the merged object
#' @export
merge.cevodata <- function(x, y, name = "Merged datasets", ...) {
  genome <- if (x$genome == y$genome) x$genome else "multiple genomes"
  metadata <- bind_rows(x$metadata, y$metadata)
  cd <- init_cevodata(name = name, genome = genome) |>
    add_sample_data(metadata)

  cd$SNVs <- bind_assays(x, y, "SNVs")
  cd$CNVs <- bind_assays(x, y, "CNVs")
  cd$clones <- bind_assays(x, y, "clones")
  cd$models <- bind_assays(x, y, "models")

  message("Setting active SNVs to ", x$active_SNV)
  message("Setting active CNVs to ", x$active_CNV)
  message("Setting active clones to ", x$active_clones)
  cd$active_SNVs <- NULL
  cd$active_CNVs <- NULL
  cd$active_clones <- NULL
  active_assays <- list(
    active_SNVs = x$active_SNVs,
    active_CNVs = x$active_CNVs,
    active_clones = x$active_clones
  )
  cd <- c(cd, active_assays)
  structure(cd, class = "cevodata")
}


bind_assays <- function(x, y, slot_name) {
  assays <- c(names(x[[slot_name]]), names(y[[slot_name]])) |>
    unique()
  if (length(assays) > 0) {
    assays |>
      set_names(assays) |>
      map(~bind_rows(x[[slot_name]][[.x]], y[[slot_name]][[.x]]))
  } else {
    list()
  }
}
