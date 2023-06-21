
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
      dentro_2015_correction() |>
      relocate("CCF", .after = "VAF")
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
      CCF = .data$u / .data$m
    )
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
#'   the CCF is found in the SNV tibble, and VAF otherwise
#' @param bins Number of interval bins
#' @param verbose Verbose?
#'
#' @return <cevodata> object
#' @export
intervalize_mutation_frequencies <- function(object,
                                             which_snvs = default_SNVs(object),
                                             column = mutation_frequencies_column_name(object, which_snvs),
                                             bins = NULL,
                                             verbose = get_cevomod_verbosity()) {
  msg("Calculating f intervals, using ", column, " column", verbose = verbose)

  snvs <- SNVs(object, which_snvs) |>
    cut_f_intervals(bins = bins, column = column)
  object |>
    add_SNV_data(snvs, name = which_snvs)
}


mutation_frequencies_column_name <- function(object, which_snvs = default_SNVs(object)) {
  snvs <- SNVs(object, which_snvs)
  if ("CCF" %in% names(snvs)) "CCF" else "VAF"
}

