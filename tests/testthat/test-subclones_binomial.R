
fit_binomial_models_cols <- c("N", "component", "cellularity", "N_mutations", "BIC")

test_that("fit_binomial_models() works with very few remaining mutations", {
  residuals_1 <- tibble(VAF = 1:100/100) |>
    mutate(
      neutral_resid_clones = case_when(
        VAF == 4/100 ~ 100,
        VAF == 6/100 ~ 50,
        TRUE ~ 0
      )
    )
  res <- fit_binomial_models(residuals_1, N = 1:3)
  expect_equal(nrow(res), 1)
  expect_named(res, fit_binomial_models_cols)
})


test_that("fit_binomial_models() works with no remaining mutations", {
  residuals_0 <- tibble(
    VAF = 1:100/100,
    neutral_resid_clones = 0
  )
  res <- fit_binomial_models(residuals_0, N = 1:3)
  expect_equal(nrow(res), 0)
  expect_named(res, fit_binomial_models_cols)
})


test_that("get_binomial_predictions() works", {
  clones <- tibble(
    component = c("Clone", "Subclone 1"),
    N_mutations = c(300, 100),
    cellularity = c(.5, .2),
    median_DP = 100
  )
  binomial <- get_binomial_predictions(clones)
  map(binomial, sum)
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
    mean_DP     = c(46, 28),
    median_DP   = c(52, 26),
    sd_DP       = c(25, 10),
  )
  # plot_clones(clones_non_overlapping)
  expect_false(any_binomial_distibutions_correlate(clones_non_overlapping))

  clones_overlapping <- tibble(
    sample_id   = "S1",
    component   = c("Clone", "Subclone 1", "Subclone 2"),
    cellularity = c(0.33, 0.17, 0.12),
    N_mutations = c(732, 866, 406),
    mean_DP     = c(54, 25, 38),
    median_DP   = c(60, 25, 37),
    sd_DP       = c(20, 8, 12),
  )
  # plot_clones(clones_overlapping)
  expect_true(any_binomial_distibutions_correlate(clones_overlapping))
})


test_that("rebinarize_distribution() does not change the number of mutations", {
  clones <- tibble(
    sample_id   = "S1",
    component   = c("Clone", "Subclone 1"),
    cellularity = c(0.33, 0.16),
    N_mutations = c(748, 1256),
    mean_DP     = c(46, 28),
    median_DP   = c(52, 26),
    sd_DP       = c(25, 10),
  )
  predictions <- get_binomial_predictions(clones)
  sums <- c(sum(predictions$Clone), sum(predictions$`Subclone 1`))
  expect_equal(sums, clones$N_mutations)
  predictions2 <- predictions |>
    select(VAF, pred = Clone) |>
    rebinarize_distribution(VAFs = 1:50/50)
  sums2 <- sum(predictions2$pred)
  expect_equal(sums2, clones$N_mutations[[1]])
})
