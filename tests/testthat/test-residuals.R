data("tcga_brca_test")


test_that("calc_powerlaw_model_residuals returns powerlaw curves with similar numbers of mutations for different numbers of bins", {
  cd <- tcga_brca_test
  cd$active_models <- "powerlaw_fixed"
  cd1 <- tcga_brca_test |>
    intervalize_mutation_frequencies() |>
    calc_SFS(bins = 50) |>
    calc_powerlaw_model_residuals("powerlaw_fixed")
  cd2 <- tcga_brca_test |>
    intervalize_mutation_frequencies() |>
    calc_SFS(bins = 100) |>
    calc_powerlaw_model_residuals("powerlaw_fixed")
  # plot_models(cd1)
  # plot_models(cd2)
  n1 <- get_residuals(cd1) |>
    filter(f > 0.1) |>
    pull(powerlaw_pred) |>
    sum()
  n2 <- get_residuals(cd2) |>
    filter(f > 0.1) |>
    pull(powerlaw_pred) |>
    sum()
  expect_true((abs(n1 - n2) / n1) < 0.02)
})


test_that("plot_sampling_rate returns correct ggplot object", {
  p <- plot_sampling_rate(tcga_brca_test)
  vdiffr::expect_doppelganger("plot-sampling-rate", p)
})


test_that("plot_residuals_powerlaw_model returns correct ggplot object", {
  p <- plot_residuals_powerlaw_model(tcga_brca_test)
  vdiffr::expect_doppelganger("plot-resuduals-powerlaw-model", p)
})


test_that("plot_residuals_full_model returns correct ggplot object", {
  p <- plot_residuals_full_model(tcga_brca_test)
  vdiffr::expect_doppelganger("plot-residuals-full-model", p)
})


test_that("plot_binomial_fits_vs_powerlaw_residuals_bars returns correct ggplot object", {
  p <- plot_binomial_fits_vs_powerlaw_residuals_bars(tcga_brca_test)
  vdiffr::expect_doppelganger("plot_binomial_fits_vs_powerlaw_residuals_bars", p)
})
