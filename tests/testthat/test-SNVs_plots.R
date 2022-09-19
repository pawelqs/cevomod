data("snvs_test")

test_that("plot_mutations() with drivers works", {
  p <- snvs_test |>
    plot_mutations(drivers = "BRCA")
  vdiffr::expect_doppelganger("plot_mutations(drivers = 'BRCA')", p)
})


test_that("plot_mutations() with gene list works", {
  p <- snvs_test |>
    plot_mutations(genes = c("TP53", "TBX3", "BRCA1"))
  vdiffr::expect_doppelganger("plot_mutations(genes = TP53,TBX3,BRCA1)", p)
})


test_that("filter_SNVs() works", {
  res <- snvs_test |>
    filter_SNVs(genes = c("TP53", "TBX3", "BRCA1"))
  expect_s3_class(res, "cevo_SNVs_tbl")
  expect_equal(nrow(res), 6)
  expect_identical(res$gene_symbol, c("TBX3", "TBX3", "BRCA1", "TP53", "BRCA1", "TBX3"))
  expect_identical(res$start, c(115109720, 115117319, 41244256, 7579882, 41244131, 115117717))
})
