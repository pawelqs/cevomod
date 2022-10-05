data("tcga_brca_test")


test_that("plot_CNV_heatmap works", {
  ht <- plot_CNV_heatmap(tcga_brca_test, "seg_mean", verbose = FALSE)
  expect_s4_class(ht, "Heatmap")
})
