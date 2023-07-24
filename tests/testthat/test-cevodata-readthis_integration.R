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



test_that("adding FACETS data works", {
  facets_dir <- system.file("extdata", "FACETS", package = "readthis")
  data <- readthis::read_facets_cnvs(facets_dir)
  cd <- init_cevodata("Test dataset") |>
    add_data(data)
  expect_s3_class(cd, "cevodata")
  expect_s3_class(CNVs(cd), "tbl")
  expect_equal(cd$active_CNVs, "FACETS")
  expect_equal(dim(CNVs(cd)), c(128, 18))
  expect_equal(cd$metadata$purity, c(0.3, 0.3))
  expect_equal(cd$metadata$purity, cd$metadata$facets_purity)
})



test_that("adding Mutect2 data works", {
  path <- system.file("extdata", "Mutect", package = "readthis")
  data <- readthis::read_mutect_snvs(
    path,
    patient_id_pattern = "(?<=Mutect\\/)[:alnum:]*(?=\\.)",
    verbose = FALSE
  )
  cd <- init_cevodata("Test dataset") |>
    add_data(data)
  expect_s3_class(cd, "cevodata")
  expect_s3_class(SNVs(cd), "tbl")
  expect_equal(cd$active_SNVs, "Mutect")
  expect_equal(dim(SNVs(cd)), c(16, 14))
  expect_equal(cd$metadata$sample_id, c("S1_L1", "S1_P1", "S2_L1", "S2_P1"))
  expect_equal(cd$metadata$patient_id, c("S1", "S1", "S2", "S2"))
})



test_that("adding Strelka data works", {
  path <- system.file("extdata", "Strelka", package = "readthis")
  data <- readthis::read_strelka_somatic_snvs(
    path,
    patient_id_pattern = "(?<=Strelka\\/)[:alnum:]*(?=\\.)",
    verbose = FALSE
  ) |>
    mutate(sample_id = str_c(patient_id, sample_id, sep = "_"))
  cd <- init_cevodata("Test dataset") |>
    add_data(data)
  expect_s3_class(cd, "cevodata")
  expect_s3_class(SNVs(cd), "tbl")
  expect_equal(cd$active_SNVs, "Strelka")
  expect_equal(dim(SNVs(cd)), c(18, 11))
  expect_equal(cd$metadata$sample_id, c("S1_TUMOR", "S2_TUMOR"))
})
