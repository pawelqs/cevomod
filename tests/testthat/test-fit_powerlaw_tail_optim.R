data("tcga_brca_test")
set_cevomod_verbosity(0)

object <- tcga_brca_test |>
  intervalize_mutation_frequencies()



test_that("fit_powerlaw_tail_optim models have non-negative objective fun value", {
  name <- "powerlaw_optim"
  bootstraps <- FALSE
  allowed_zero_bins <- 2
  y_treshold <- 1
  y_threshold_pct <- 0.01
  av_filter <- c(1/3, 1/3, 1/3)
  peak_detection_upper_limit <- 0.3
  reward_upper_limit <- 0.4
  control <- list(maxit = 1000, ndeps = c(0.1, 0.01))
  verbose <- get_cevomod_verbosity()

  object <- fit_powerlaw_tail_optim(object)
  td <- get_models(object, "powerlaw_optim")
  expect_true(all(td$value > 1990))
  expect_s3_class(get_powerlaw_models(object, "powerlaw_optim"), "cevo_powerlaw_models")
})



test_that("fit_powerlaw_tail_optim() bootstrapping test", {
  name <- "powerlaw_optim"
  bootstraps <- 2
  allowed_zero_bins <- 2
  y_treshold <- 1
  y_threshold_pct <- 0.01
  av_filter <- c(1/3, 1/3, 1/3)
  peak_detection_upper_limit <- 0.3
  reward_upper_limit <- 0.4
  control <- list(maxit = 1000, ndeps = c(0.1, 0.01))
  verbose <- TRUE
  object <- object |>
    filter(sample_id %in% object$metadata$sample_id[1:2])

  suppressWarnings({
    object <- fit_powerlaw_tail_optim(object, bootstraps = bootstraps)
  })

  bootstrap_models <- get_models(object, "powerlaw_optim_bootstraps")
  expect_equal(nrow(bootstrap_models), 4)
  expect_true(all(bootstrap_models$value > 1990))
  expect_s3_class(bootstrap_models, "cevo_powerlaw_models")
  expect_s3_class(bootstrap_models, "cevo_bootstrap_powerlaw_models")

  models <- get_models(object, "powerlaw_optim")
  expected_columns <- c(
    "sample_id", "model", "component",
    "A", "A.lower", "A.upper", "alpha", "alpha.lower", "alpha.upper"
  )
  expect_named(models, expected_columns)
  expect_equal(nrow(models), 2)
  expect_s3_class(models, "cevo_powerlaw_models")
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

