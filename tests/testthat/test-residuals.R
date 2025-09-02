data("tcga_brca_fitted")


test_that("Fitting powerlaw_fixed model returns the same residuals", {
  coefs <- test_path("testdata", "tcga_brca_2samples.coefs_powerlaw_fixed.tsv") |>
    read_tsv(show_col_types = FALSE) |>
    filter(best)
  sfs <- test_path("testdata", "tcga_brca_2samples.SFS.tsv") |>
    read_tsv(show_col_types = FALSE)
  expected <- test_path("testdata", "tcga_brca_2samples.residuals_powerlaw_fixed.tsv") |>
    read_tsv(show_col_types = FALSE)

  residuals <- calc_powerlaw_model_residuals(coefs, sfs)
  expect_equal(residuals, expected)
})



test_that("calc_powerlaw_model_residuals returns powerlaw curves with similar numbers of mutations for different numbers of bins", {
  cd <- tcga_brca_fitted |>
    intervalize_mutation_frequencies()
  active_models(cd) <- "powerlaw_fixed"

  resid1 <- calc_powerlaw_model_residuals(
    get_model_coefficients(cd),
    cd |>
      calc_SFS(bins = 50) |>
      get_SFS()
  )
  resid2 <- calc_powerlaw_model_residuals(
    get_model_coefficients(cd),
    cd |>
      calc_SFS(bins = 100) |>
      get_SFS()
  )
  n1 <- resid1 |>
    filter(f > 0.1) |>
    pull(powerlaw_pred) |>
    sum()
  n2 <- resid2 |>
    filter(f > 0.1) |>
    pull(powerlaw_pred) |>
    sum()
  expect_true((abs(n1 - n2) / n1) < 0.02)
})


test_that("plot_sampling_rate returns correct ggplot object", {
  p <- plot_sampling_rate(tcga_brca_fitted)
  vdiffr::expect_doppelganger("plot-sampling-rate", p)
})


test_that("plot_residuals_powerlaw_model returns correct ggplot object", {
  p <- plot_residuals_powerlaw_model(tcga_brca_fitted)
  vdiffr::expect_doppelganger("plot-resuduals-powerlaw-model", p)
})


test_that("plot_residuals_full_model returns correct ggplot object", {
  p <- plot_residuals_full_model(tcga_brca_fitted)
  vdiffr::expect_doppelganger("plot-residuals-full-model", p)
})


test_that("plot_binomial_fits_vs_powerlaw_residuals_bars returns correct ggplot object", {
  p <- plot_binomial_fits_vs_powerlaw_residuals_bars(tcga_brca_fitted)
  vdiffr::expect_doppelganger("plot_binomial_fits_vs_powerlaw_residuals_bars", p)
})
