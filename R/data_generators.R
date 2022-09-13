
#' Generate neutral SNVs
#'
#' @param mut_rate a in f(x) = a/x^2 formula
#' @param sample_below mutations below this value will be sampled
#' @return SNVs tibble
#' @export
#'
#' @examples
#' generate_neutral_snvs()
generate_neutral_snvs <- function(mut_rate = 2, sample_below = 0.15) {
  tibble(
    sample_id = "S1",
    VAF = seq(.01, 1, by = .01),
    n = floor(mut_rate / .data$VAF^2)
  ) |>
    mutate(
      n = if_else(
        .data$VAF <= sample_below,
        floor(mut_rate / sample_below^2 * (.data$VAF / sample_below)^2),
        n
      ),
      mut_id = map(.data$n, ~ tibble(mut_id = rep("a", times = .x)))
    ) |>
    unnest(.data$mut_id)
}
