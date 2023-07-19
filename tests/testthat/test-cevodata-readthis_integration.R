test_that("adding ASCAT data works", {
  ascat_dir <- system.file("extdata", "ASCAT", package = "readthis")
  data <- readthis::read_ascat_files(ascat_dir, sample_id_pattern = "(?<=ASCAT\\/)[:alnum:]*(?=\\.)")
  cd <- init_cevodata("Test dataset") |>
    add_data(data)
  expect_s3_class(cd, "cevodata")
  expect_s3_class(CNVs(cd), "tbl")
  expect_equal(cd$active_CNVs, "ASCAT")
  expect_equal(dim(CNVs(cd)), c(20, 8))
  expect_equal(cd$metadata$purity, c(0.99322, 0.99322))
  expect_equal(cd$metadata$purity, cd$metadata$ascat_purity)
})
