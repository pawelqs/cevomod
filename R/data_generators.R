
#' Generate neutral SNVs
#'
#' @param mut_rate a in f(x) = a/x^2 formula
#' @param sample_below mutations below this value will be sampled
#' @param resolution resolution
#' @param DP sequncing depth
#' @return SNVs tibble
#' @export
#'
#' @examples
#' generate_neutral_snvs()
generate_neutral_snvs <- function(mut_rate = 2, sample_below = 0.15, resolution = 0.01, DP = 100) {
  tibble(
    patient_id = "S1",
    sample_id = "S1",
    VAF = seq(.01, 1, by = resolution),
    DP = DP,
    n = floor(mut_rate / .data$VAF^2)
  ) |>
    mutate(
      n = if_else(
        .data$VAF <= sample_below,
        floor(mut_rate / sample_below^2 * (.data$VAF / sample_below)^2),
        n
      ),
      alt_reads = round(.data$DP * .data$VAF),
      ref_reads = round(.data$DP * (1 - .data$VAF)),
      mutation_id = map(.data$n, ~ tibble(mutation_id = rep("a", times = .x)))
    ) |>
    unnest("mutation_id") |>
    mutate(
      mutation_id = str_c("chr1", row_number())
    )
}
