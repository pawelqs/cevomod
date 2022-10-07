data("tcga_brca_test")
snvs <- SNVs(tcga_brca_test)
cnvs <- CNVs(tcga_brca_test)
meta <- tcga_brca_test$metadata


test_that("print.cevodata runs without error works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_sample_data(meta)
  output <- print(cd) |>
    capture.output()
  expect_true(length(cd) > 1)
})


test_that("print.cevodata runs without error foe empty object", {
  cd <- init_cevodata("TCGA BRCA small")
  output <- print(cd) |>
    capture.output()
  expect_true(length(cd) > 1)
})
