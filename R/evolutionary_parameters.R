

#' Get evolutionary parameters from the model
#'
#' Code of those functions is only re-formatted code from MOBSTER R package
#' by Caravagna, Williams et al. https://github.com/caravagnalab/mobster
#'
#' @param object cevodata object or models tibble
#' @param ... other arguments
#' @param models_name models_name
#' @name evo_params
NULL


## ------------------------ Mutation rates by Williams ------------------------

#' @describeIn evo_params Get mutation rates by Williams
#' Use ...
#' @export
get_mutation_rates <- function(object, ...) {
  UseMethod("get_mutation_rates")
}


#' @describeIn evo_params Get mutation rates by Williams
#' @export
get_mutation_rates.cevodata <- function(object, models_name = "williams_neutral", ...) {
  mutation_rates <- get_models(object, models_name) |>
    filter(.data$component == "Neutral tail") |>
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
  #   summarise(N = sum(powerlaw_pred)) |>
  #   left_join(bin_widths) |>
  #   mutate(
  #     min_VAF = min_VAF - bin_width / 2,
  #     max_VAF = max_VAF + bin_width / 2,
  #     mu = N / ( 1/(2*min_VAF) - 1/(2*max_VAF) ) / 2
    # )
}


#' @describeIn evo_params Get mutation rates by Williams
#' @export
get_mutation_rates.tbl_df <- function(object, ...) {
  require_columns(object, "sample_id", "component", "A")
  mutation_rates <- object |>
    filter(.data$component == "Neutral tail") |>
    transmute(
      .data$sample_id,
      mutation_rate_williams = .data$A
    )
  mutation_rates
}


## --------------------- Selection coefs by Williams --------------------------
# Code below is re-formatted code from MOBSTER R package
# by Caravagna, Williams et al. https://github.com/caravagnalab/mobster


#' @describeIn evo_params Get subclonal selection coefficients Williams
#'
#' Use properties of subclone fit to calculate selection intensity, selection
#' is defined as the relative growth rates of host tumour cell
#' populations (\eqn{\lambda h}) vs subclone (\eqn{\lambda s}):
#' \deqn{1+s=\lambda h / \lambda s}
#'
#' @param Nmax Time when tumour is sampled (in tumour doublings)
#' @export
get_selection_coefficients <- function(object, ...) {
  UseMethod("get_selection_coefficients")
}


#' @describeIn evo_params Get subclonal selection coefficients Williams
#' @export
get_selection_coefficients.cevodata <- function(
          object,
          models_name = "williams_neutral_subclones",
          Nmax = 10^10, ...) {
  get_models(object, models_name) |>
    get_selection_coefficients()
}


#' @describeIn evo_params Get subclonal selection coefficients Williams
#' @export
get_selection_coefficients.tbl_df <- function(object, Nmax = 10^10, ...) {
  require_columns(object, "sample_id", "component", "A", "N_mutations", "cellularity")
  mutation_rates <- get_mutation_rates(object)

  subclones <- object |>
    filter(str_detect(.data$component, "Subclone")) |>
    drop_na_columns() |>
    select("sample_id", "component", "N_mutations", "cellularity") |>
    mutate(subclone_frequency = 2 * .data$cellularity) |> # need ccf so times by 2
    select(-"cellularity")

  dt <- subclones |>
    left_join(mutation_rates, by = "sample_id") |>
    mutate(
      emergence_time = get_emergence_time(.data$N_mutations, .data$mutation_rate_williams)
    ) |>
    nest(subclones = c("component", "N_mutations", "subclone_frequency", "emergence_time"))

  dt |>
    rowwise("sample_id") |>
    reframe(mobster_evolutionary_parameters(.data$subclones, .data$mu)) |>
    left_join(mutation_rates, by = "sample_id")
}


# Emergence time is negative, when N is not big enough, eg. comparable with mu
get_emergence_time <- function(N, mu) {
  (N / mu) / (2 * log(2)) - (-digamma(1) / log(2))
}


mobster_evolutionary_parameters <- function(subclones, mu, Nmax = 10^10) {
  nsubclones <- nrow(subclones)

  subclones$selection <- if (nsubclones == 1) {
    time_end <- log(Nmax * (1 - subclones$subclone_frequency)) / log(2)
    selection(subclones$emergence_time, time_end, subclones$subclone_frequency)
  } else if (nsubclones == 2) {
    if (are_subclones_nested(subclones)) { # pigeon hole principle
      largestsubclone <- max(subclones$subclone_frequency)
      time_end <- log(Nmax * (1 - largestsubclone)) / log(2)
      selection2clonenested(subclones$emergence_time, time_end, subclones$subclone_frequency)
    }
    else {
      time_end <- log(Nmax * ( 1 - subclones$subclone_frequency[1] - subclones$subclone_frequency[2] )) / log(2)
      selection2clone(subclones$emergence_time, time_end, subclones$subclone_frequency)
    }
  }
  return(subclones)
}


are_subclones_nested <- function(subclones) {
  sum(subclones$subclone_frequency) > 1
}


# Calculate strength of selection of subclone
# @param time time subclone emerges (in tumour doublings)
# @param time_end Time when tumour is sampled (in tumour doublings)
# @param subclonefrequency Frequency of subclones
selection <- function(time, time_end, subclonefrequency) {
  x1 <- log(2) * time
  x2 <- log(subclonefrequency / (1 - subclonefrequency))
  x3 <- log(2) * (time_end - time)

  s <- ((x1 + x2) / x3)
  return(s)
}


# Calculate strength of selection for 2 independent subclones
# @param times emergence times (in tumour doublings)
# @param time_end Time when tumour is sampled (in tumour doublings)
# @param frequencies Frequencies of both subclones
selection2clone <- function(times, time_end, frequencies) {
  time1 <- times[1]
  time2 <- times[2]
  subclonefrequency1 <- frequencies[1]
  subclonefrequency2 <- frequencies[2]

  x1 <- log(2) * time1
  x2 <-
    log(subclonefrequency1 / (1 - subclonefrequency1 - subclonefrequency2))
  x3 <- log(2) * (time_end - time1)
  s1 <- ((x1 + x2) / x3)

  x1 <- log(2) * time2
  x2 <-
    log(subclonefrequency2 / (1 - subclonefrequency1 - subclonefrequency2))
  x3 <- log(2) * (time_end - time2)
  s2 <- ((x1 + x2) / x3)

  return(c(s1, s2))
}


# Calculate strength of selection for 2 nested subclones
selection2clonenested <- function(times, time_end, frequencies) {
  time1 <- times[1]
  time2 <- times[2]
  subclonefrequency1 <- frequencies[1]
  subclonefrequency2 <- frequencies[2]

  x1 <- log(2) * time1
  x2 <- log((subclonefrequency1 - subclonefrequency2)
            / (1 - subclonefrequency1))
  x3 <- log(2) * (time_end - time1) # I replaced time with time1, I am sure it was a typo; PK
  s1 <- ((x1 + x2) / x3)

  x1 <- log(2) * time2
  x2 <- log((subclonefrequency2)
            / (1 - subclonefrequency1))
  x3 <- log(2) * (time_end - time2)
  s2 <- ((x1 + x2) / x3)

  return(c(s1, s2))
}

