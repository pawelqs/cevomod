data("tcga_brca_test")
snvs <- SNVs(tcga_brca_test)
cnvs <- CNVs(tcga_brca_test)
meta <- tcga_brca_test$metadata


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


test_that("Adding SNVs to cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_SNV_data(snvs = snvs, "head")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$SNVs), 2)
  expect_s3_class(cd$SNVs[[1]], "tbl_df")
  expect_s3_class(cd$SNVs[[2]], "tbl_df")
  expect_equal(cd$active_SNVs, "head")
})


test_that("Adding SNVs to cevodata extends metadata", {
  small_snvs <- snvs |>
    filter(sample_id %in% c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"))
  cd <- init_cevodata("TCGA BRCA small", snvs = small_snvs)
  expect_equal(cd$metadata$sample_id, c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"))
  cd <- add_SNV_data(cd, snvs = snvs, "head")
  expect_equal(nrow(cd$metadata), 4)
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
  expect_equal(SNVs(cd, "snvs") |> nrow(), 21570)
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


test_that("Adding CNVs to cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", snvs = snvs) |>
    add_CNV_data(cnvs, "tcga")
  expect_s3_class(cd, "cevodata")
  expect_equal(length(cd$CNVs), 1)
  expect_s3_class(cd$CNVs[[1]], "tbl_df")
  expect_equal(cd$active_CNVs, "tcga")
})


test_that("Adding CNVs to cevodata extends metadata", {
  small_cnvs <- cnvs |>
    filter(sample_id %in% c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"))
  cd <- init_cevodata("TCGA BRCA small", cnvs = small_cnvs)
  expect_equal(cd$metadata$sample_id, c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"))
  cd <- add_CNV_data(cd, cnvs = cnvs, "head")
  expect_equal(nrow(cd$metadata), 4)
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
  expect_equal(CNVs(cd, "cnvs") |> nrow(), 714)
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


test_that("Setting active SNV assay on cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnvs = cnvs) |>
    add_CNV_data(cnvs = cnvs, "tcga2")
  default_CNVs(cd) <- "cnvs"
  expect_equal(default_CNVs(cd), "cnvs")
  expect_error(default_CNVs(cd) <- "xxx")
})


test_that("Adding sample data to cevodata works", {
  meta <- tibble(
    sample_id =  c("TCGA-AC-A23H-01", "TCGA-AN-A046-01"),
    meta1 = 1:2
  )
  cd <- init_cevodata("TCGA BRCA small") |>
    add_sample_data(meta)
  expect_equal(cd$metadata$sample_id, meta$sample_id)
  expect_equal(cd$metadata$meta1, meta$meta1)

  meta2 <- tibble(
    sample_id =  c("TCGA-D8-A27V-01", "TCGA-AN-A046-01"),
    meta2 = c("a", "b")
  )
  cd <- add_sample_data(cd, meta2)
  expect_equal(cd$metadata$sample_id, c("TCGA-AC-A23H-01", "TCGA-AN-A046-01", "TCGA-D8-A27V-01"))
  expect_equal(cd$metadata$meta1, c(1, 2, NA_real_))
  expect_equal(cd$metadata$meta2, c(NA_character_, "b", "a"))

  cd <- add_SNV_data(cd, snvs)
  expect_equal(nrow(cd$metadata), 4)
  expect_named(cd$metadata, c("sample_id", "meta1", "meta2"))
  expect_true(all(unique(snvs$sample_id) %in% cd$metadata$sample_id))
  expect_equal(cd$metadata$meta1, c(1, 2, NA_real_, NA_real_))
})
