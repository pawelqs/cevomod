data("tcga_brca_test")


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
