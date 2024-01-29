## ----------------------- Mutation counts by component -----------------------

#' Get numbers of neutral, clonal, subclonal ect variants
#' @param object cevodata object
#' @param models_name models name
#' @param include_filtered include the component of the filtered mutations
#' @export
count_mutations_by_component <- function(object, models_name = active_models(object), include_filtered = FALSE) {
  neutral_tail_counts <- count_neutral_tail_mutations(object, models_name)
  mut_counts <- object |>
    get_models(models_name) |>
    select("sample_id", "component", "N_mutations") |>
    mutate(component = str_replace(.data$component, "powerlaw", "Neutral tail")) |>
    left_join(neutral_tail_counts, by = c("sample_id", "component")) |>
    mutate(
      N_mutations = if_else(.data$component == "Neutral tail", .data$N, .data$N_mutations)
    )

  if (include_filtered) {
    filtered_mut_counts <- count_neutral_tail_filtered_mutations(object, models_name) |>
      rename(N_mutations = "N") |>
      filter(.data$sample_id %in% mut_counts$sample_id)
    mut_counts <- bind_rows(mut_counts, filtered_mut_counts) |>
      arrange(.data$sample_id)
  }

  mut_counts <- mut_counts |>
    mutate(N = sum(.data$N_mutations), .by = "sample_id") |>
    mutate(frac = .data$N_mutations / .data$N)
  mut_counts
}


#' @describeIn count_mutations_by_component Count neutral tail mutations
#' @export
count_neutral_tail_mutations <- function(object, models_name = active_models(object)) {
  object |>
    get_residuals(models_name = models_name) |>
    mutate(neutral_tail_muts = pmin(.data$SFS, .data$powerlaw_pred)) |>
    summarise(
      component = "Neutral tail",
      N = sum(.data$neutral_tail_muts) |> round() |> as.integer(),
      .by = "sample_id"
    ) |>
    complete(
      sample_id = object$metadata$sample_id,
      fill = lst(component = "Neutral tail")
    )
}


#' @describeIn count_mutations_by_component Count filtered mutations with VAF higher than 0.01
#' @export
count_neutral_tail_filtered_mutations <- function(object, models_name = active_models(object)) {
  object |>
    calc_SFS(bins = 100) |>
    calc_powerlaw_model_residuals(models_name) |>
    get_residuals(models_name) |>
    filter(.data$VAF > 0.01) |>
    mutate(
      filtered_muts = if_else(.data$powerlaw_pred > .data$SFS, .data$powerlaw_pred - .data$SFS, 0)
    ) |>
    summarise(
      component = "Filtered mutations",
      N = sum(.data$filtered_muts) |> round() |> as.integer(),
      .by = "sample_id"
    ) |>
    complete(
      sample_id = object$metadata$sample_id,
      fill = lst(component = "Filtered mutations")
    )
}
