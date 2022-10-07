data("tcga_brca_test")


test_that("plot_sampling_rate returns correct ggplot object", {
  p <- plot_sampling_rate(tcga_brca_test)
  vdiffr::expect_doppelganger("plot-sampling-rate", p)
})


test_that("plot_residuals_neutral_model returns correct ggplot object", {
  p <- plot_residuals_neutral_model(tcga_brca_test)
  vdiffr::expect_doppelganger("plot-resuduals-neutral-model", p)
})


test_that("plot_residuals_full_model returns correct ggplot object", {
  p <- plot_residuals_full_model(tcga_brca_test)
  vdiffr::expect_doppelganger("plot-residuals-full-model", p)
})
