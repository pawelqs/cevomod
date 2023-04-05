data("tcga_brca_test")

test_that("plot_mutations() with drivers works", {
  p <- tcga_brca_test |>
    plot_mutations(drivers = "BRCA")
  vdiffr::expect_doppelganger("plot_mutations(drivers = 'BRCA')", p)
})


test_that("plot_mutations() with gene list works", {
  p <- tcga_brca_test |>
    plot_mutations(genes = c("TP53", "TBX3", "BRCA1"))
  p2 <- tcga_brca_test |>
    SNVs() |>
    plot_mutations(genes = c("TP53", "TBX3", "BRCA1"))
  vdiffr::expect_doppelganger("plot_mutations(genes = TP53,TBX3,BRCA1)", p)
  vdiffr::expect_doppelganger("plot_mutations(genes = TP53,TBX3,BRCA1)", p2)
})


test_that("filter_SNVs() works", {
  res <- SNVs(tcga_brca_test) |>
    filter_SNVs(genes = c("TP53", "TBX3", "BRCA1"))
  expect_s3_class(res, "cevo_snvs")
  expect_equal(nrow(res), 5)
  expect_identical(res$gene_symbol, c("TBX3", "TBX3", "BRCA1", "TP53", "TBX3"))
  expect_identical(res$pos, c(115109720, 115117319, 41244256, 7579882, 115117717))
})
