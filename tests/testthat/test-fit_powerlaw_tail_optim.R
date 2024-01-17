data("tcga_brca_fitted")
verbose::verbose(cevoverse = 0)


# --------------------------- fit_powerlaw_tail_optim -------------------------

test_that("fit_powerlaw_tail_optim() returns models with non-negative objective values", {
  name <- "powerlaw_optim"
  bootstraps <- FALSE
  allowed_zero_bins <- 2
  y_treshold <- 1
  y_threshold_pct <- 0.01
  av_filter <- c(1/3, 1/3, 1/3)
  peak_detection_upper_limit <- 0.3
  reward_upper_limit <- 0.4
  control <- list(maxit = 1000, ndeps = c(0.1, 0.01))
  verbose <- get_verbosity()

  object <- tcga_brca_fitted |>
    intervalize_mutation_frequencies() |>
    fit_powerlaw_tail_optim()
  td <- get_model_coefficients(object, "powerlaw_optim")
  expect_true(all(td$value > 1990))
  expect_s3_class(get_models(object, "powerlaw_optim"), "cevo_powerlaw_models")
})


# test_that("fit_powerlaw_tail_optim() bootstrapping test", {
#   name <- "powerlaw_optim"
#   bootstraps <- 2
#   allowed_zero_bins <- 2
#   y_treshold <- 1
#   y_threshold_pct <- 0.01
#   av_filter <- c(1/3, 1/3, 1/3)
#   peak_detection_upper_limit <- 0.3
#   reward_upper_limit <- 0.4
#   control <- list(maxit = 1000, ndeps = c(0.1, 0.01))
#   verbose <- TRUE
#   object <- object |>
#     filter(sample_id %in% object$metadata$sample_id[1:2])
#
#   suppressWarnings({
#     object <- fit_powerlaw_tail_optim(object, bootstraps = bootstraps)
#   })
#
#   bootstrap_models <- get_models(object, "powerlaw_optim_bootstraps")
#   expect_equal(nrow(bootstrap_models), 4)
#   expect_true(all(bootstrap_models$value > 1990))
#   expect_s3_class(bootstrap_models, "cevo_powerlaw_models")
#   expect_s3_class(bootstrap_models, "cevo_bootstrap_powerlaw_models")
#
#   models <- get_models(object, "powerlaw_optim")
#   expected_columns <- c(
#     "sample_id", "model", "component",
#     "A", "A.lower", "A.upper", "alpha", "alpha.lower", "alpha.upper"
#   )
#   expect_named(models, expected_columns)
#   expect_equal(nrow(models), 2)
#   expect_s3_class(models, "cevo_powerlaw_models")
# })


# --------------------------------- Helpers ------------------------------------

test_that("td_objective_function() works", {
  object <- tcga_brca_fitted
  sfs <- get_SFS(object)
  bounds <- get_f_range(SNVs(object), pct_left = 0.02, pct_right = 0.98)
  nbins <- SNVs(object) |>
    group_by(.data$sample_id) |>
    summarise(nbins = median(DP))

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


test_that("get_non_zero_SFS_range() works", {
  SFS <- tibble(
    sample_id = "S1",
    f = 1:100,
    y = 10
  ) |>
    mutate(
      y = case_when(
        f < 12 ~ 0,
        f == 40 ~ 0,
        f %in% c(54, 55) ~ 1,
        f > 65 ~ 0,
        TRUE ~ y
      ),
      f = f / 100
    )
  allowed_zero_bins <- 1
  y_treshold <- 1

  expected <- tibble(
    sample_id = "S1",
    from = 0.12,
    to = 0.65
  )
  res <- get_non_zero_SFS_range(SFS)

  expect_identical(res, expected)
})
