set_cevomod_verbosity(0)

object <- test_data |>
  prepare_SNVs() |>
  fit_powerlaw_tail_optim()

N <- 1:3
powerlaw_model_name <- "powerlaw_optim"
upper_f_limit <- 0.75
snvs_name <- "snvs"


test_that("fit_subclones() mclust works", {
  res <- object |>
    fit_subclones(method = "mclust")
  expect_false(
    is.null(res$models$powerlaw_optim_subclones)
  )
})


test_that("fit_subclones() bmix works", {
  suppressWarnings(
    suppressMessages({
      res <- object |>
        fit_subclones(method = "BMix")
    })
  )
  expect_false(
    is.null(res$models$powerlaw_optim_subclones)
  )
})


fit_binomial_models_cols <- c("N", "component", "cellularity", "N_mutations", "BIC")


test_that("fit_binomial_models() works with very few remaining mutations", {
  residuals_1 <- tibble(f = 1:100/100) |>
    mutate(
      powerlaw_resid_clones = case_when(
        f == 4/100 ~ 100,
        f == 6/100 ~ 50,
        TRUE ~ 0
      )
    )
  res <- fit_binomial_models(residuals_1, N = 1:3)
  expect_equal(nrow(res), 1)
  expect_named(res, fit_binomial_models_cols)
})


test_that("fit_binomial_models() works with no remaining mutations", {
  residuals_0 <- tibble(
    f = 1:100/100,
    powerlaw_resid_clones = 0
  )
  res <- fit_binomial_models(residuals_0, N = 1:3)
  expect_equal(nrow(res), 0)
  expect_named(res, fit_binomial_models_cols)
})


test_that("get_binomial_predictions() works", {
  intervals <- tibble(
    f = 1:100/100,
    f_interval = cut_number(f, n = 100)
  )
  clones <- tibble(
    component = c("Clone", "Subclone 1"),
    N_mutations = c(300, 100),
    cellularity = c(.5, .2),
    sequencing_DP = 100
  )
  binomial <- get_binomial_predictions(clones, intervals)
  expect_equal(sum(binomial$Clone), 300)
  expect_equal(sum(binomial$`Subclone 1`), 100)
  expect_equal(max(binomial$Clone), 23.876771, tolerance = 0.0001)
  expect_equal(max(binomial$`Subclone 1`), 9.930021, tolerance = 0.0001)
})


test_that("any_binomial_distibutions_correlate() works", {
  clones_non_overlapping <- tibble(
    sample_id   = "S1",
    component   = c("Clone", "Subclone 1"),
    cellularity = c(0.33, 0.16),
    N_mutations = c(748, 1256),
    sequencing_DP   = c(52, 26)
  )
  # plot_clones(clones_non_overlapping)
  expect_false(any_binomial_distibutions_correlate(clones_non_overlapping))

  clones_overlapping <- tibble(
    sample_id   = "S1",
    component   = c("Clone", "Subclone 1", "Subclone 2"),
    cellularity = c(0.33, 0.17, 0.12),
    N_mutations = c(732, 866, 406),
    sequencing_DP   = c(60, 25, 37)
  )
  # plot_clones(clones_overlapping)
  expect_true(any_binomial_distibutions_correlate(clones_overlapping))
})


test_that("rebinarize_distribution() does not change the number of mutations", {
  predictions <- tibble(
    f = 1:10/10,
    Clone = 1:10,
    `Subclone 1` = c(0, 3, 5, 20, 15, 13, 9, 0, 0, 0)
  )
  predictions2 <- rebinarize_distribution(predictions, n_bins = 25)
  expect_equal(sum(predictions$Clone), sum(predictions2$Clone))
  expect_equal(sum(predictions$`Subclone 1`), sum(predictions2$`Subclone 1`))
  expect_equal(nrow(predictions2), 25)
})


