data("tcga_brca_test")


test_that("plot_CNV_heatmap works", {
  ht <- plot_CNV_heatmap(tcga_brca_test, "seg_mean", verbose = FALSE)
  ht <- ComplexHeatmap::draw(ht)
  vdiffr::expect_doppelganger("plot_cnv_heatmap_tcga_brca_test", ht)
})
