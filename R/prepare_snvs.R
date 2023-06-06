
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


## ----------------------------------- MCF -------------------------------------


#' Calc mutation frequencies
#'
#' Currently no method is implemented and this function only initializes
#' f column in active SNVs slot with values from VAF column
#'
#' @name mutation_frequencies


#' @rdname mutation_frequencies
#' @export
calc_mutation_frequencies <- function(object, ...) {
  UseMethod("calc_mutation_frequencies")
}


#' @rdname mutation_frequencies
#' @param object <cevodata> object
#' @param method Available methods:
#'   - use_VAF - use raw VAF avalues as mutation frequencies
#'   - Dentro - use CNV-based frequency correction from
#'   [Dentro et al. 'Principles of Reconstructing the Subclonal Architecture of Cancers' (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5538405/)
#'   doi: 10.1101/cshperspect.a026625
#' @param which_snvs Which SNVs to use
#' @param which_cnvs Which CNVs to use
#' @param rm_intermediate_cols Should the columns used to get MCF be removed?
#' @param ... Other arguments (not used)
#' @export
calc_mutation_frequencies.cevodata <- function(object,
                                               method = "Dentro",
                                               which_snvs = default_SNVs(object),
                                               which_cnvs = default_CNVs(object),
                                               rm_intermediate_cols = TRUE,
                                               ...) {
  if (method == "use_VAF") {
    snvs <- SNVs(object, which_snvs) |>
      mutate(MCF = .data$VAF, .after = "VAF")
    intermediate_cols <- c()
  } else if (method == "Dentro") {
    cnvs <- CNVs(object, which_cnvs) |>
      select("sample_id", "chrom", "start", "end", "total_cn", "normal_cn")
    purities <- get_purities(object)
    snvs <- SNVs(object, which_snvs) |>
      join_CNVs(cnvs) |>
      left_join(purities, by = "sample_id") |>
      dentro_2015_correction() |>
      relocate("MCF", .after = "VAF")
    intermediate_cols <- c("start", "end", "total_cn", "normal_cn", "purity", "u", "m")
  } else {
    stop("Currently supported methods: 'use_VAF' and 'Dentro'")
  }

  if (rm_intermediate_cols) {
    snvs <- snvs |>
      select(-any_of(intermediate_cols))
  }

  object |>
    add_SNV_data(snvs, name = which_snvs)
}


#' @describeIn mutation_frequencies CNV-based frequency correction method described by Dentro et al. (2015)
#'
#' Method was described in [Dentro et al. 'Principles of Reconstructing the Subclonal Architecture of Cancers' (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5538405/)
#' doi: 10.1101/cshperspect.a026625
#'
#' @param tbl tibble that contains columns: VAF, total_cn, normal_cn, purity
#' @return The same tibble with new columns: u, m, f. f contains the corrected
#'   Mutated Cell Frequencies
#' @export
dentro_2015_correction <- function(tbl) {
  tbl |>
    mutate(
      u = .data$VAF * (1 / .data$purity) * (.data$purity * .data$total_cn + (1 - .data$purity) * .data$normal_cn),
      m = if_else(.data$u < 1, 1, round(.data$u)),
      MCF = .data$u / .data$m
    )
}


## --------------------------------- Cut intervals -----------------------------


#' Intervalize mutation frequencies
#' @param object cevodata obj
#' @param which_snvs which SNVs to use
#' @param bins number of VAF interval bins
#' @export
intervalize_mutation_frequencies <- function(object,
                                             which_snvs = default_SNVs(object),
                                             bins = NULL) {
  snvs <- SNVs(object, which_snvs) |>
    cut_f_intervals(bins = bins)
  object |>
    add_SNV_data(snvs, name = which_snvs)
}
