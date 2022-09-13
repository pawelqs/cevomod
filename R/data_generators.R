
generate_neutral_snvs <- function(mut_rate = 2, sample_below = 0.15) {
  tibble(
    sample_id = "S1",
    VAF = seq(.01, 1, by = .01),
    n = floor(mut_rate/VAF^2)
  ) |>
  mutate(
    n = if_else(VAF <= sample_below, floor(mut_rate/sample_below^2 * (VAF/sample_below)^2), n),
    mut_id = map(n, ~tibble(mut_id = rep("a", times = .x)))
  ) |>
  unnest(mut_id)
}
