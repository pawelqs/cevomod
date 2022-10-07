data("tcga_brca_test")
snvs <- SNVs(tcga_brca_test)
cnvs <- CNVs(tcga_brca_test)
meta <- tcga_brca_test$metadata


test_that("Filtering cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnvs = cnvs) |>
    add_CNV_data(cnvs = cnvs, "tcga2") |>
    add_SNV_data(snvs = snvs, "tcga") |>
    filter(sample_id == "TCGA-AC-A23H-01")
  expect_equal(nrow(cd$CNVs$cnvs), 285)
  expect_equal(nrow(cd$CNVs$tcga2), 285)
  expect_equal(nrow(cd$SNVs$tcga), 6419)
  expect_equal(unique(cd$CNVs$cnvs$sample_id), "TCGA-AC-A23H-01")
  expect_equal(unique(cd$CNVs$tcga2$sample_id), "TCGA-AC-A23H-01")
  expect_equal(unique(cd$SNVs$tcga$sample_id), "TCGA-AC-A23H-01")
})


test_that("Merging cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnvs = cnvs) |>
    add_CNV_data(cnvs = cnvs, "tcga2") |>
    add_SNV_data(snvs = snvs, "tcga")
  cd1 <- filter(cd, sample_id == "TCGA-AC-A23H-01")
  cd2 <- filter(cd, sample_id != "TCGA-AC-A23H-01")
  cd_merged <- merge(cd1, cd2, name = "TCGA BRCA small", verbose = FALSE)
  output <- print(cd_merged) |>
    capture.output()
  expect_equal(object.size(cd_merged), object.size(cd))
})
