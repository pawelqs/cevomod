

#' Get mutation rates by Williams
#' @param object cevodata object  with fitted cevomod models
#' @export
get_mutation_rates <- function(object) {
  mutation_rates <- get_neutral_models(object) |>
    transmute(
      .data$sample_id,
      mutation_rate_williams = .data$A
    )
  mutation_rates
  # models <- get_models(object)
  #
  # residuals <- get_residuals(object)
  # bin_widths <- residuals |>
  #   group_by(sample_id) |>
  #   summarise(bin_width = get_interval_width(VAF_interval))
  #
  # min_VAF <- 0.2
  # max_VAF <- 0.8
  # residuals |>
  #   filter(VAF > min_VAF, VAF < max_VAF) |>
  #   group_by(sample_id) |>
  #   summarise(N = sum(neutral_pred)) |>
  #   left_join(bin_widths) |>
  #   mutate(
  #     min_VAF = min_VAF - bin_width / 2,
  #     max_VAF = max_VAF + bin_width / 2,
  #     mu = N / ( 1/(2*min_VAF) - 1/(2*max_VAF) ) / 2
    # )
}


#' Get subclonal selection coefficients Williams
#' @param object cevodata object
#' @param Nmax Time when tumour is sampled (in tumour doublings)
#' @export
get_selection_coefficients <- function(object, Nmax = 10^10) {
  mutation_rates <- get_mutation_rates(object)

  subclones <- get_models(object) |>
    filter(str_detect(component, "Subclone")) |>
    drop_na_columns() |>
    select("sample_id", "component", "N_mutations", "cellularity") |>
    mutate(cellularity = 2 * .data$cellularity)

  evolutionary_parameters <- subclones |>
    left_join(mutation_rates, by = "sample_id") |>
    group_by(.data$sample_id) |>
    mutate(
      subclones_in_sample = str_detect(.data$component, "Subclone") |> sum(),
      emergence_time = get_emergence_time(.data$N_mutations, .data$mutation_rate_williams),
      sequencing_time = get_sequencing_time(.data$subclones_in_sample, Nmax, .data$cellularity)
    ) |>
    rowwise() |>
    mutate(
      selection = calc_selection_coef(.data$emergence_time, .data$sequencing_time, .data$cellularity)
    ) |>
    ungroup() |>
    select("sample_id", "component", "emergence_time", "sequencing_time", "selection")

  evolutionary_parameters
}


get_emergence_time <- function(N, mu) {
  (N / mu) / (2 * log(2)) - (-digamma(1) / log(2))
}


get_sequencing_time <- function(n_subclones, N_max, cellularity) {
  if_else(
    n_subclones == 1,
    log(N_max * (1 - cellularity)) / log(2),
    NA_real_
  )
}


# Calculate strength of selection of subclone (from MOBSTER)
#
# Use properties of subclone fit to calculate selection intensity, selection
# is defined as the relative growth rates of host tumour cell
# populations (\eqn{\lambda h}) vs subclone (\eqn{\lambda s}):
# \deqn{1+s=\lambda h / \lambda s}
#
# @param time time subclone emerges (in tumour doublings)
# @param time_end Time when tumour is sampled (in tumour doublings)
# @param subclonefrequency Frequency of subclones
calc_selection_coef <- function(time, time_end, subclonefrequency) {
  x1 <- log(2) * time
  x2 <- log(subclonefrequency / (1 - subclonefrequency))
  x3 <- log(2) * (time_end - time)
  s <- ((x1 + x2) / x3)
  s
}
