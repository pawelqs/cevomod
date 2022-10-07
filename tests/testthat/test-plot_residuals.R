data("tcga_brca_test")


test_that("plot_sampling_rate returns correct ggplot object", {
  p <- plot_sampling_rate(tcga_brca_test)
  expect_s3_class(p, "ggplot")
})


test_that("plot_residuals_neutral_model returns correct ggplot object", {
  p <- plot_residuals_neutral_model(tcga_brca_test)
  expect_s3_class(p, "ggplot")
})


test_that("plot_residuals_full_model returns correct ggplot object", {
  p <- plot_residuals_full_model(tcga_brca_test)
  expect_s3_class(p, "ggplot")
})
