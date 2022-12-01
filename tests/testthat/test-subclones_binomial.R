
cd <- tcga_brca_test |>
  filter(sample_id == "TCGA-AN-A046-01")

test_that("predict_binoms works", {
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

  residuals_0 <- tibble(
    VAF = 1:100/100,
    neutral_resid_clones = 0
  )
  res <- fit_binomial_models(residuals_0, N = 1:3)
  expect_equal(nrow(res), 0)
})


clones_non_overlapping <- tibble(
  sample_id   = "S1",
  component   = c("Clone", "Subclone 1"),
  cellularity = c(0.33, 0.16),
  N_mutations = c(748, 1256),
  mean_DP     = c(46, 28),
  median_DP   = c(52, 26),
  sd_DP       = c(25, 10),
)

clones_overlapping <- tibble(
  sample_id   = "S1",
  component   = c("Clone", "Subclone 1", "Subclone 2"),
  cellularity = c(0.33, 0.17, 0.12),
  N_mutations = c(732, 866, 406),
  mean_DP     = c(54, 25, 38),
  median_DP   = c(60, 25, 37),
  sd_DP       = c(20, 8, 12),
)

# clones_overlapping |>
#   rename(sequencing_DP = .data$median_DP) |>
#   pmap(get_binomial_distribution) |>
#   bind_rows(.id = "component") |>
#   ggplot(aes(VAF, pred)) +
#   geom_point()
# plot(x$`Subclone 2`, x$Clone)


test_that("any_binomial_distibutions_correlate() works", {
  expect_false(any_binomial_distibutions_correlate(clones_non_overlapping))
  expect_true(any_binomial_distibutions_correlate(clones_overlapping))
})

