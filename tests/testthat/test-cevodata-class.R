data("tcga_brca_test")
snvs <- SNVs(tcga_brca_test)
cnvs <- CNVs(tcga_brca_test)


test_that("init_cevodata works", {
  cd <- init_cevodata("TCGA BRCA small")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$SNVs), 0)
  expect_equal(length(cd$cNVs), 0)
  expect_equal(cd$genome, "unknown")
  cd <- init_cevodata("TCGA BRCA small", genome = "hg19")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$SNVs), 0)
  expect_equal(length(cd$cNVs), 0)
  expect_equal(cd$genome, "hg19")
})


test_that("init_cevodata with SNV assay works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs)
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$SNVs), 1)
  expect_s3_class(cd$SNVs[[1]], "tbl_df")
  expect_equal(cd$active_SNVs, "snvs")
})


test_that("print.cevodata runs without error works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs)
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


test_that("Adding SNV to cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_SNV_data(snvs = snvs, "head")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$SNVs), 2)
  expect_s3_class(cd$SNVs[[1]], "tbl_df")
  expect_s3_class(cd$SNVs[[2]], "tbl_df")
  expect_equal(cd$active_SNVs, "head")
})


test_that("Adding incorrect SNV to cevodata throws error", {
  expect_error(init_cevodata("TCGA BRCA small", snvs = cnvs))
})


test_that("Getting SNV from cevodata works", {
  tcga2 <- head(snvs)
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_SNV_data(snvs = tcga2, "head")
  expect_identical(SNVs(cd), tcga2)
  expect_equal(SNVs(cd) |> nrow(), 6)
  expect_identical(SNVs(cd, "head"), tcga2)
  expect_equal(SNVs(cd, "snvs") |> nrow(), 23492)
  expect_error(SNVs(cd, "xxx"))
})


test_that("Getting active SNV assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_SNV_data(snvs = snvs, "head")
  expect_equal(default_SNVs(cd), "head")
})


test_that("Setting active SNV assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_SNV_data(snvs = snvs, "head")
  default_SNVs(cd) <- "snvs"
  expect_equal(default_SNVs(cd), "snvs")
  expect_error(default_SNVs(cd) <- "xxx")
})


test_that("Adding CNV to cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_CNV_data(cnvs, "tcga")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$CNVs), 1)
  expect_s3_class(cd$CNVs[[1]], "tbl_df")
  expect_equal(cd$active_CNVs, "tcga")
})


test_that("Adding incorrect CNV to cevodata throws error", {
  expect_error(init_cevodata("TCGA BRCA small", cnvs = snvs))
})


test_that("Getting CNV from cevodata works", {
  cnvs2 <- head(cnvs)
  cd <- init_cevodata("TCGA BRCA small", cnvs = cnvs) |>
    add_CNV_data(cnvs2, "tcga2")
  expect_identical(CNVs(cd), cnvs2)
  expect_equal(CNVs(cd) |> nrow(), 6)
  expect_identical(CNVs(cd, "tcga2"), cnvs2)
  expect_identical(CNVs(cd, "cnvs"), cnvs)
  expect_equal(CNVs(cd, "cnvs") |> nrow(), 762)
  expect_error(CNVs(cd, "xxx"))
})


test_that("Getting active CNV assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnvs = cnvs) |>
    add_CNV_data(cnvs = cnvs, "tcga")
  expect_equal(default_CNVs(cd), "tcga")
})


test_that("Setting active SNV assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnvs = cnvs) |>
    add_CNV_data(cnvs = cnvs, "tcga2")
  default_CNVs(cd) <- "cnvs"
  expect_equal(default_CNVs(cd), "cnvs")
  expect_error(default_CNVs(cd) <- "xxx")
})


test_that("Filtering cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnvs = cnvs) |>
    add_CNV_data(cnvs = cnvs, "tcga2") |>
    add_SNV_data(snvs = snvs, "tcga") |>
    filter(patient_id == "TCGA-AC-A23H-01")
  expect_equal(nrow(cd$CNVs$cnvs), 285)
  expect_equal(nrow(cd$CNVs$tcga2), 285)
  expect_equal(nrow(cd$SNVs$tcga), 6419)
  expect_equal(unique(cd$CNVs$cnvs$patient_id), "TCGA-AC-A23H-01")
  expect_equal(unique(cd$CNVs$tcga2$patient_id), "TCGA-AC-A23H-01")
  expect_equal(unique(cd$SNVs$tcga$patient_id), "TCGA-AC-A23H-01")
})
