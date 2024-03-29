
## --------------------------------- Prepare SNVs ------------------------------

#' Prepare SNVs for analyses
#' @param object cevodata obj
#' @param bins number of VAF interval bins
#' @param verbose verbose?
#' @export
prepare_SNVs <- function(object,
                         bins = NULL,
                         verbose = TRUE) {
  object |>
    calc_mutation_frequencies(method = "Dentro") |>
    intervalize_mutation_frequencies(bins = bins) |>
    calc_SFS()
}


## ----------------------------------- CCF -------------------------------------


#' Calc mutation frequencies
#'
#' Calculates the CNV-corrected mutation frequencies. Implemented methods:
#' - Dentro - calculates the Cancer Cell Fraction (CCF) using the formulas from
#'   [Dentro et al. *Principles of Reconstructing the Subclonal Architecture of Cancers* (2015)](https://doi.org/10.1101/cshperspect.a026625)
#'
#' @param object <cevodata> object
#' @param method Available methods: Dentro
#' @param which_snvs Which SNVs to use
#' @param which_cnvs Which CNVs to use
#' @param rm_intermediate_cols Should the columns used to get CCF be removed?
#' @param verbose Verbose?
#'
#' @return <cevodata> object
#' @export
calc_mutation_frequencies <- function(object,
                                      method = "Dentro",
                                      which_snvs = default_SNVs(object),
                                      which_cnvs = default_CNVs(object),
                                      rm_intermediate_cols = TRUE,
                                      verbose = get_cevomod_verbosity()) {
  if (method == "Dentro") {
    cnvs <- CNVs(object, which_cnvs) |>
      select("sample_id", "chrom", "start", "end", "total_cn", "normal_cn")
    purities <- get_purities(object)
    snvs <- SNVs(object, which_snvs) |>
      join_CNVs(cnvs) |>
      left_join(purities, by = "sample_id") |>
      dentro_2015_correction()
    intermediate_cols <- c("start", "end", "total_cn", "normal_cn", "purity", "u", "m")
  } else {
    stop("Currently supported methods: 'Dentro'")
  }

  if (rm_intermediate_cols) {
    snvs <- snvs |>
      select(-any_of(intermediate_cols))
  }

  nas <- snvs |>
    filter(is.na(.data$CCF)) |>
    nrow()
  nas_pct <- round(nas*100/nrow(snvs), digits = 2)
  msg(nas, " variants (", nas_pct, " %), have NA CCF value", verbose = verbose)

  object |>
    add_SNV_data(snvs, name = which_snvs)
}


#' @describeIn calc_mutation_frequencies Implements the CNV-based frequency
#' correction method described in [Dentro et al. 'Principles of Reconstructing the Subclonal Architecture of Cancers' (2015)](https://doi.org/10.1101/cshperspect.a026625)
#'
#' @param tbl tibble that contains columns: VAF, total_cn, normal_cn, purity
#'
#' @return The same tibble with new columns: u, m, CCF (Cancer Cell Fraction)
#' @export
dentro_2015_correction <- function(tbl) {
  tbl |>
    mutate(
      u = .data$VAF * (1 / .data$purity) * (.data$purity * .data$total_cn + (1 - .data$purity) * .data$normal_cn),
      m = if_else(.data$u < 1, 1, round(.data$u)),
      CCF = .data$u / .data$m,
      `CCF/2` = .data$CCF / 2
    ) |>
    relocate("CCF", "CCF/2", .after = "VAF")
}


## --------------------------------- Cut intervals -----------------------------


#' Intervalize the mutation frequencies
#'
#' Intervalize the mutation frequencies for the subsequent analyses and plots.
#' Adds f_interval column to the SNV tibble
#'
#' @param object <cevodata> object
#' @param which_snvs Which SNVs to use
#' @param column Which frequency measure column to intervalize? By default, uses
#'   the CCF/2  if it is found in the SNV tibble, and VAF otherwise
#' @param bins Number of interval bins
#' @param ... Other args
#' @param verbose Verbose?
#'
#' @return <cevodata> object
#' @export
intervalize_mutation_frequencies <- function(object, ...) {
  UseMethod("intervalize_mutation_frequencies")
}



#' @rdname intervalize_mutation_frequencies
#' @export
intervalize_mutation_frequencies.cevodata <- function(object,
                                                      which_snvs = default_SNVs(object),
                                                      column = get_frequency_measure_name(object, which_snvs),
                                                      bins = NULL,
                                                      verbose = get_cevomod_verbosity(),
                                                      ...) {
  snvs <- SNVs(object, which_snvs) |>
    intervalize_mutation_frequencies(column, bins, verbose)
  object |>
    add_SNV_data(snvs, name = which_snvs)
}


#' @rdname intervalize_mutation_frequencies
#' @export
intervalize_mutation_frequencies.cevo_snvs <- function(object,
                                                       column = get_frequency_measure_name(object),
                                                       bins = NULL,
                                                       verbose = get_cevomod_verbosity(),
                                                       ...) {
  msg("Calculating f intervals, using ", column, " column", verbose = verbose)
  object |>
    cut_f_intervals(bins = bins, column = column) |>
    as_cevo_snvs()
}


#' Decide which mutation frequency measure to use
#'
#' Used during the intervalization of mutation frequencies. Uses CCF/2 if found
#' in the object, and VAF otherwise
#'
#' @param object object
#' @param which_snvs Which SNVs to use?
#' @param ... other params (not used now)
#'
#' @export
get_frequency_measure_name <- function(object, ...) {
  UseMethod("get_frequency_measure_name")
}


#' @rdname get_frequency_measure_name
#' @export
get_frequency_measure_name.cevodata <- function(object, which_snvs = default_SNVs(object), ...) {
  SNVs(object, which_snvs) |>
    get_frequency_measure_name()
}


#' @rdname get_frequency_measure_name
#' @export
get_frequency_measure_name.cevo_snvs <- function(object, ...) {
  if ("CCF/2" %in% names(object)) "CCF/2" else "VAF"
}
