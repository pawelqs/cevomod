## ----------------------- Mutation counts by component -----------------------

#' Get numbers of neutral, clonal, subclonal ect variants
#' @param object cevodata object
#' @param models_name models name
#' @param include_missing include the component of the missing mutations
#' @param min_f f threshold - power law tail goes to the Infinity, so we need
#'   to set a reasonable threshold
#' @export
count_mutations_by_component <- function(object, models_name = active_models(object), include_missing = FALSE, min_f = 0.01) {
  neutral_tail_counts <- count_powerlaw_tail_mutations(object, models_name)
  mut_counts <- object |>
    get_model_coefficients(models_name) |>
    select("sample_id", "component", "N_mutations") |>
    mutate(component = str_replace(.data$component, "powerlaw tail", "Powerlaw tail")) |>
    left_join(neutral_tail_counts, by = c("sample_id", "component")) |>
    mutate(
      N_mutations = if_else(.data$component == "Powerlaw tail", .data$N, .data$N_mutations)
    ) |>
    select(-"N")

  if (include_missing) {
    filtered_mut_counts <- count_missing_powerlaw_tail_mutations(object, models_name, min_f = min_f)
    if (!is.null(filtered_mut_counts)) {
      filtered_mut_counts <- filtered_mut_counts |>
        rename(N_mutations = "N") |>
        filter(.data$sample_id %in% mut_counts$sample_id)
      mut_counts <- bind_rows(mut_counts, filtered_mut_counts) |>
        arrange(.data$sample_id)
    }
  }

  mut_counts <- mut_counts |>
    mutate(N = sum(.data$N_mutations), .by = "sample_id") |>
    mutate(frac = .data$N_mutations / .data$N)
  mut_counts
}


#' @describeIn count_mutations_by_component Count neutral tail mutations
#' @export
count_powerlaw_tail_mutations <- function(object, models_name = active_models(object)) {
  object |>
    get_model_residuals(model_name = models_name) |>
    mutate(powerlaw_tail_muts = pmin(.data$SFS, .data$powerlaw_pred)) |>
    summarise(
      component = "Powerlaw tail",
      N = sum(.data$powerlaw_tail_muts) |> round() |> as.integer(),
      .by = "sample_id"
    ) |>
    complete(
      sample_id = get_metadata(object)$sample_id,
      fill = lst(component = "Powerlaw tail")
    )
}


#' @describeIn count_mutations_by_component Count filtered mutations with f higher than threshold
#' @export
count_missing_powerlaw_tail_mutations <- function(object,
                                                  models_name = active_models(object),
                                                  min_f = 0.01) {
  residuals <- get_model_residuals(object, models_name)

  bin_numbers <- count(residuals, .data$sample_id)
  bin_numbers_equal <- all(bin_numbers$n == bin_numbers$n[1])
  if (!bin_numbers_equal) {
    warning(
      "The number of bins is not equal for all samples.",
      "Numbers of missing powerlaw tail mutations will not be estimated."
    )
    return(NULL)
  }

  residuals |>
    filter(.data$f > min_f) |>
    mutate(
      missing_mutations = if_else(.data$powerlaw_pred > .data$SFS, .data$powerlaw_pred - .data$SFS, 0)
    ) |>
    summarise(
      component = "Missing mutations",
      N = sum(.data$missing_mutations) |> round() |> as.integer(),
      .by = "sample_id"
    ) |>
    complete(
      sample_id = get_metadata(object)$sample_id,
      fill = lst(component = "Missing mutations")
    )
}
