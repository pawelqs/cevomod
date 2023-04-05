
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
  patient_ids <- unique(new_object$metadata[["patient_id"]])

  new_object$SNVs   <- map(new_object$SNVs,   ~filter(.x, sample_id %in% ids))
  new_object$CNVs   <- map(new_object$CNVs,   ~filter(.x, sample_id %in% ids))
  new_object$models <- map(new_object$models, ~filter(.x, sample_id %in% ids))
  new_object$misc   <- map(new_object$misc,   ~filter(.x, sample_id %in% ids))
  new_object$misc_by_sample  <- map(new_object$misc_by_sample,  ~.x[ids])
  new_object$misc_by_patient <- map(new_object$misc_by_patient, ~.x[patient_ids])

  if (is_cevodata_singlepatient(new_object)) {
    class(new_object) <- c("singlepatient_cevodata", class(new_object))
  }
  new_object
}


filter_joined_models <- function(joined_models, ...) {
  if (is.null(joined_models)) {
    return(NULL)
  }
  samples <- joined_models |>
    map(~tibble(sample = c(.x$rowsample, .x$colsample))) |>
    bind_rows(.id = "patient_id")

  samples_kept <- samples |>
    filter(...) |>
    mutate(sample_kept = TRUE)

  kept_patients <- samples |>
    left_join(samples_kept, by = c("patient_id", "sample")) |>
    group_by(.data$patient_id) |>
    filter(all(.data$sample_kept)) |>
    pull("patient_id") |>
    unique()

  joined_models[kept_patients]
}


#' Merge two cevodata objects
#' @inheritParams base::merge
#' @param name Name of the merged object
#' @param verbose Show messages?
#' @param .id datasets names will be saved to this metadata column, if provided
#' @export
merge.cevodata <- function(x, y, name = "Merged datasets", verbose = TRUE, .id = NULL, ...) {
  genome <- if (x$genome == y$genome) x$genome else "multiple genomes"
  if (!is.null(.id)) {
    x$metadata[[.id]] <- x$name
    y$metadata[[.id]] <- y$name
  }
  metadata <- bind_rows(x$metadata, y$metadata)
  cd <- init_cevodata(name = name, genome = genome) |>
    add_sample_data(metadata)

  cd$SNVs <- bind_assays(x, y, "SNVs")
  cd$CNVs <- bind_assays(x, y, "CNVs")
  cd$models <- bind_assays(x, y, "models")
  cd$misc <- bind_assays(x, y, "misc")
  cd$misc_by_sample <- map2(x$misc_by_sample, y$misc_by_sample, ~union(.x, .y))
  cd$misc_by_patient <- map2(x$misc_by_patient, y$misc_by_patient, ~union(.x, .y))

  if (verbose) {
    message("Setting active SNVs to ", x$active_SNV)
    message("Setting active CNVs to ", x$active_CNV)
  }
  cd$active_SNVs <- NULL
  cd$active_CNVs <- NULL
  active_assays <- list(
    active_SNVs = x$active_SNVs,
    active_CNVs = x$active_CNVs
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


#' Split object
#' @param object object to split
#' @param ... other arguments
#' @export
split_by <- function(object, ...) {
  UseMethod("split_by")
}


#' @describeIn split_by Split cevodata object
#' @param object cevodata object
#' @param x name of column in metadata
#' @export
split_by.cevodata <- function(object, x, ...) {
  split_names <- object$metadata[[x]] |>
    unique()
  splits <- split_names |>
    set_names(split_names) |>
    map(~filter(object, .data[[x]] == .x))
  class(splits) <- c("cevo_splits", "list")
  splits
}


# @export
# print.cevo_splits <- function(x, ...) {
#   cli::cat_line("<cevo_splits> object. Splits:")
#   cli::cat_line(paste0(names(x), collapse = ", "))
# }


#' Update cevodata object with values from another object
#' @param object object to update
#' @param object2 object to use
#' @param ... other args, unused now
#' @export
update.cevodata <- function(object, object2, ...) {
  sample_ids <- union(object$metadata$sample_id, object2$metadata$sample_id)
  object <- object |>
    filter(.data$sample_id %not in% object2$metadata$sample_id) |>
    merge(object2)
  object$metadata <- object$metadata |>
    mutate(sample_id = parse_factor(.data$sample_id, levels = sample_ids)) |>
    arrange(.data$sample_id) |>
    mutate(across("sample_id", as.character))
  object
}
