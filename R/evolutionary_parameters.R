
## ---------------------- Get evolutionary parameters --------------------------

#' Get evolutionary parameters from the model
#'
#' Code of those functions is only re-formatted code from MOBSTER R package
#' by Caravagna, Williams et al. https://github.com/caravagnalab/mobster
#'
#'
#'
#' Use properties of subclone fit to calculate selection intensity, selection
#' is defined as the relative growth rates of host tumour cell
#' populations (\eqn{\lambda h}) vs subclone (\eqn{\lambda s}):
#' \deqn{1+s=\lambda h / \lambda s}
#'
#' @param Nmax Time when tumour is sampled (in tumour doublings)
#'
#' @param object cevodata object or models tibble
#' @param ... other arguments
#' @param models_name models_name
#' @name evo_params
NULL


#' @describeIn evo_params Get evolutionary parameters
#' @params object cevodata or models object
#' @export
get_evolutionary_parameters <- function(object, ...) {
  UseMethod("get_evolutionary_parameters")
}


#' @export
get_evolutionary_parameters.cevodata <- function(
    object,
    models_name = active_models(object),
    Nmax = 10^10,
    ...) {
  object |>
    get_models(models_name) |>
    get_evolutionary_parameters(Nmax = Nmax)
}


#' @export
get_evolutionary_parameters.cv_powerlaw_models <- function(object, ...) {
  get_mutation_rates(object$coefs)
}


#' @export
get_evolutionary_parameters.cv_powerlaw_subclones_models <- function(
    object,
    Nmax = 10^10,
    ...) {
  models <- object
  mutation_rates <- get_mutation_rates(models$coefs)
  subclones <- models$coefs |>
    filter(str_detect(.data$component, "Subclone")) |>
    select("sample_id", "component", "N_mutations", "frequency") |>
    mutate(cellular_frequency = 2 * .data$frequency) |> # need ccf so times by 2
    select(-"frequency")

  dt <- subclones |>
    left_join(mutation_rates, by = "sample_id") |>
    mutate(
      emergence_time = mobster_emergence_time(.data$N_mutations, .data$mutation_rate)
    ) |>
    nest(subclones = c("component", "N_mutations", "cellular_frequency", "emergence_time"))

  dt |>
    rowwise("sample_id", "mutation_rate") |>
    reframe(mobster_evolutionary_params(.data$subclones, Nmax = Nmax))
}


## ------------------------------ Mutation rates -------------------------------

get_mutation_rates <- function(coefs) {
  require_columns(coefs, "sample_id", "component", "A")
  mutation_rates <- coefs |>
    filter(.data$component %in% c("Neutral tail", "powerlaw tail")) |>
    transmute(
      .data$sample_id,
      mutation_rate = .data$A,
      mutation_rate_method = "williams"
    )
  mutation_rates
}


# get_mutation_rates_like_MOBSTER <- function() {
  # models <- get_models(object)
  # residuals <- get_residuals(object)
  # bin_widths <- residuals |>
  #   group_by(sample_id) |>
  #   summarise(bin_width = get_interval_width(f_interval))
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
# }



## --------------------- Selection coefs by Williams --------------------------
# Code below is re-formatted code from MOBSTER R package
# by Caravagna, Williams et al. https://github.com/caravagnalab/mobster


mobster_emergence_time <- function(N, mu) {
  # Emergence time is negative, when N is not big enough, eg. comparable with mu
  (N / mu) / (2 * log(2))# - (-digamma(1) / log(2))
}


mobster_evolutionary_params <- function(subclones, Nmax = 10^10) {
  require_columns(subclones, "component", "N_mutations", "cellular_frequency", "emergence_time")
  nsubclones <- nrow(subclones)

  if (nsubclones == 1) {
    time_end <- log(Nmax * (1 - subclones$cellular_frequency)) / log(2)
    s <- selection(subclones$emergence_time, time_end, subclones$cellular_frequency)
  } else if (nsubclones == 2) {
    if (are_subclones_nested(subclones)) { # pigeon hole principle
      largestsubclone <- max(subclones$cellular_frequency)
      time_end <- log(Nmax * (1 - largestsubclone)) / log(2)
      s <- selection2clonenested(subclones$emergence_time, time_end, subclones$cellular_frequency)
    }
    else {
      time_end <- log(Nmax * ( 1 - subclones$cellular_frequency[1] - subclones$cellular_frequency[2] )) / log(2)
      s <- selection2clone(subclones$emergence_time, time_end, subclones$cellular_frequency)
    }
  }

  subclones$time_end <- time_end
  subclones$selection_coef <- s
  return(subclones)
}


are_subclones_nested <- function(subclones) {
  sum(subclones$cellular_frequency) > 1
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
  x2 <- log(subclonefrequency1 / (1 - subclonefrequency1 - subclonefrequency2))
  x3 <- log(2) * (time_end - time1)
  s1 <- ((x1 + x2) / x3)

  x1 <- log(2) * time2
  x2 <- log(subclonefrequency2 / (1 - subclonefrequency1 - subclonefrequency2))
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
  x2 <- log((subclonefrequency1 - subclonefrequency2) / (1 - subclonefrequency1))
  x3 <- log(2) * (time_end - time1) # I replaced time with time1, sure it was a typo; PK
  s1 <- ((x1 + x2) / x3)

  x1 <- log(2) * time2
  x2 <- log((subclonefrequency2) / (1 - subclonefrequency1))
  x3 <- log(2) * (time_end - time2)
  s2 <- ((x1 + x2) / x3)

  return(c(s1, s2))
}
