data("tcga_brca_test")
set_cevomod_verbosity(0)

object <- tcga_brca_test |>
  intervalize_mutation_frequencies()


test_that("fit_powerlaw_tail_optim models have non-negative objective fun value", {
  object <- fit_powerlaw_tail_optim(object)
  td <- get_models(object, "powerlaw_optim")
  expect_true(all(td$value > 1990))
  expect_s3_class(get_powerlaw_models(object, "powerlaw_optim"), "cevo_powerlaw_models")
})


test_that("Testing td_objective_function", {
  object <- tcga_brca_test
  sfs <- get_SFS(object)
  bounds <- get_f_range(SNVs(object), pct_left = 0.02, pct_right = 0.98)
  nbins <- get_sample_sequencing_depths(SNVs(object)) |>
    transmute(.data$sample_id, nbins = .data$median_DP)

  data <- sfs |>
    left_join(bounds, by = "sample_id") |>
    filter(.data$f > .data$lower_bound, .data$f < .data$higher_bound) |>
    select("sample_id", "f", "y") |>
    nest_by(.data$sample_id) |>
    left_join(nbins, by = "sample_id")

  dt <- data$data[[1]]
  x <- dt$f
  y <- dt$y
  params <- c(255.9222929 / 74, 2.912159)
  res <- -td_objective_function(params, x, y)
  expect_true(res > 2000)
})
