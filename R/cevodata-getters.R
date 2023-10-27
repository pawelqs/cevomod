

## ------------------------------------ SNVs ----------------------------------

#' @rdname assays
#' @export
SNVs <- function(object, ...) {
  UseMethod("SNVs")
}

#' @return tibble
#' @describeIn assays Get SNVs from cevodata dataset
#' @export
SNVs.cevodata <- function(object, which = object$active_SNVs, ...) {
  if (which %not in% names(object$SNVs)) {
    stop(str_c(which, " does not exist in object$SNVs"))
  }
  object$SNVs[[which]]
}


#' @describeIn assays Get SNVs in the wide table form
#' @param object cevodata object
#' @param fill_na fill missing with this value
#' @export
get_SNVs_wider <- function(object, fill_na = NULL) {
  patients_to_samples <- object$metadata |>
    select("patient_id":"sample")

  snvs <- SNVs(object) |>
    select("sample_id", "chrom":"alt", "VAF") |>
    left_join(patients_to_samples, by = "sample_id") |>
    unite_mutation_id() |>
    select(-"sample_id") |>
    select("patient_id", everything()) |>
    pivot_wider(names_from = "sample", values_from = "VAF")

  if (!is.null(fill_na)) {
    snvs[is.na(snvs)] <- fill_na
  }
  snvs
}


get_SNVs_2d_matrix <- function(object,
                               rows_sample = NULL, cols_sample = NULL,
                               bins = NULL, verbose = TRUE) {
  patients_to_samples <- object$metadata |>
    select("patient_id":"sample")

  if (is.null(rows_sample) || is.null(cols_sample)) {
    rows_sample <- patients_to_samples$sample[[1]]
    cols_sample <- patients_to_samples$sample[[2]]
    msg("Using '", rows_sample, "' as rows and '", cols_sample, "' as cols", verbose = verbose)
  }
  if (n_distinct(patients_to_samples$patient_id) > 1) {
    stop("This function works only for single sample objects")
  }
  if (rows_sample %not in% patients_to_samples$sample) {
    stop("Sample requested for rows is not present for this patient")
  }
  if (cols_sample %not in% patients_to_samples$sample) {
    stop("Sample requested for cols is not present for this patient")
  }

  breaks <- object |>
    SNVs() |>
    get_interval_breaks(bins = bins)
  sample_ids <- patients_to_samples$sample_id |>
    set_names(patients_to_samples$sample)
  rowsample_breaks <- breaks[[sample_ids[rows_sample]]]
  colsample_breaks <- breaks[[sample_ids[cols_sample]]]

  # Prepare SNVs wider tibble
  mutations <-get_SNVs_wider(object, fill_na = 0)
  mutations <- mutations[c("mutation_id", rows_sample, cols_sample)]
  colnames(mutations) <- c("mutation_id", "rows_sample", "cols_sample")

  # Cut intervals
  mutations <- mutations |>
    mutate(
      rows_sample = cut(rows_sample, breaks = rowsample_breaks),
      cols_sample = cut(cols_sample, breaks = colsample_breaks)
    )
  row_intervals <- levels(mutations$rows_sample)
  col_intervals <- levels(mutations$cols_sample)

  # Prepare square matrix
  incomplete_mat <- mutations |>
    mutate(across(c("rows_sample", "cols_sample"), as.character)) |>
    group_by(.data$rows_sample, .data$cols_sample) |>
    count() |>
    ungroup() |>
    pivot_wider(names_from = "cols_sample", values_from = "n") |>
    column_to_rownames("rows_sample")
  incomplete_mat[is.na(incomplete_mat)] <- 0

  mat <- matrix(0, nrow = length(row_intervals), ncol = length(col_intervals))
  rownames(mat) <- row_intervals
  colnames(mat) <- col_intervals
  mat[rownames(incomplete_mat), colnames(incomplete_mat)] <- as.matrix(incomplete_mat)

  attr(mat, "rows_sample") <- rows_sample
  attr(mat, "cols_sample") <- cols_sample
  attr(mat, "rows_sample_id") <- sample_ids[[rows_sample]]
  attr(mat, "cols_sample_id") <- sample_ids[[cols_sample]]
  mat
}


get_SNVs_wider_intervals <- function(object, fill_na = NULL, bins = NULL) {
  metadata <- object$metadata
  if (n_distinct(metadata$patient_id) > 1) {
    stop("This function works only for single sample objects")
  }

  breaks <- object |>
    SNVs() |>
    get_interval_breaks(bins = bins)
  breaks <- breaks[metadata$sample_id]
  names(breaks) <- metadata$sample

  mutations <-get_SNVs_wider(object, fill_na = 0)
  sample_intervals <- list()
  for (sample in metadata$sample) {
    VAF_intervals <- cut(mutations[[sample]], breaks = breaks[[sample]])
    sample_intervals[[sample]] <- levels(VAF_intervals)
    mutations[[sample]] <- as.character(VAF_intervals)
  }

  attr(mutations, "sample_intervals") <- sample_intervals
  mutations
}


#' Get SNVs with merged CNVs
#' @param object cevodata object with SNVs and CNVs
#' @export
SNVs_CNVs <- function(object) {
  SNVs(object) |>
    join_CNVs(CNVs(object))
}


join_CNVs <- function(snvs, cnvs) {
  left_join(
    snvs, cnvs,
    by = join_by("sample_id", "chrom", "pos" >= "start", "pos" <= "end"),
    relationship = "many-to-one"
  )
}


## ------------------------------- Models ------------------------------------


#' Get model names
#' @param object object
#' @export
get_model_names <- function(object) {
  UseMethod("get_model_names")
}


#' @export
get_model_names.cevodata <- function(object) {
  names(object$models)
}


#' Get models from the object
#' @param object object to get the models from
#' @param ... other arguments
#' @export
get_models <- function(object, ...) {
  UseMethod("get_models")
}


#' @describeIn get_models Get models from cevodata object
#' @param which `chr` which models to get
#' @param best_only `lgl` return only the best models?
#' @export
get_models.cevodata <- function(object,
                                which = active_models(object),
                                best_only = TRUE,
                                ...) {
  models <- object$models[[which]]
  models_class <- class(models)
  if (is.null(models)) {
    stop("Slot ", which, " is empty! Fit apropriate model first")
  }
  if (best_only && !is.null(models[["best"]])) {
    models <- filter(models, .data$best)
    class(models) <- models_class
    models
  } else {
    models
  }
}


get_powerlaw_models <- function(object,
                                which = active_models(object),
                                best_only = TRUE,
                                ...) {
  models <- get_models(object, which)
  if ("cevo_powerlaw_models" %not in% class(models)) {
    stop(
      which, " is not a powerlaw model, required to fit subclones.",
      "Use another model"
    )
  }
  models
}


#' Get the name of the active models
#' @param object cebodata object
#' @export
active_models <- function(object) {
  if (is.null(object$active_models) | length(object$models) == 0) {
    stop("No models has been fitted yet!")
  }
  object$active_models
}


#' Fill N_mutations column for powerlaw models
#' @param models tibble from get_models()
#' @param cd cevodata object
#' @param models_name models name
#' @export
fix_powerlaw_N_mutations <- function(models, cd, models_name) {
  mut_counts <- count_mutations_by_component(cd, models_name, include_filtered = TRUE)
  areas_under_curves <- mut_counts |>
    filter(.data$component %in% c("Neutral tail", "Filtered mutations")) |>
    summarise(
      component = "powerlaw",
      area_under_curve = sum(.data$N_mutations),
      .by = "sample_id"
    )
  models |>
    left_join(areas_under_curves, by = c("sample_id", "component")) |>
    mutate(N_mutations = if_else(.data$component == "powerlaw", .data$area_under_curve, .data$N_mutations))
}


## ---------------------------------- CNVs -----------------------------------

#' @rdname assays
#' @export
CNVs <- function(object, ...) {
  UseMethod("CNVs")
}

#' @describeIn assays Get CNVs from cevodata dataset
#' @export
CNVs.cevodata <- function(object, which = object$active_CNVs, ...) {
  if (which %not in% names(object$CNVs)) {
    stop(str_c(which, " does not exist in object$CNVs"))
  }
  object$CNVs[[which]]
}


#' Get names of CNV variables
#' @param object object
#' @param ... other arguments
#' @export
get_CNVs_var_names <- function(object, ...) {
  UseMethod("get_CNVs_var_names")
}

#' @describeIn get_CNVs_var_names Get CNV variable names from cevodata object
#' @param which CNV assay to use
#' @export
get_CNVs_var_names.cevodata <- function(object, which = default_CNVs(object), ...) {
  cnvs_metadata <- CNVs(object, which = which) |>
    select(-"sample_id", -"chrom", -"start", -"end") |>
    drop_na_columns()
  colnames(cnvs_metadata)
}


## ---------------------------------- Other -----------------------------------

get_purities <- function(cd) {
  cd$metadata |>
    select("sample_id", "purity")
}


get_patients_data <- function(metadata) {
  patient_data_cols <- metadata |>
    group_by(.data$patient_id) |>
    summarise_all(n_distinct) |>
    map(~all(.x == 1)) |>
    keep(~.x) |>
    names()
  metadata |>
    select("patient_id", all_of(patient_data_cols))
}


get_sample_ids <- function(cd) {
  cd$metadata$sample_id
}


get_patient_sex <- function(cd) {
  cd$metadata |>
    select("sample_id", "sex")
}


#' Get sample metadata
#' @param cd cevodata object
#' @export
get_metadata <- function(cd) {
  cd$metadata
}
