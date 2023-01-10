test_that("Data generator generates correct data", {
  snvs <- generate_neutral_snvs()
  expect_s3_class(snvs, "tbl_df")
  expect_equal(nrow(snvs), 1532)
  expected_columns <- c(
    "patient_id", "sample_id", "sample", "chrom", "pos", "gene_symbol", "ref",
    "alt", "ref_reads", "alt_reads", "impact", "VAF", "DP", "n", "mut_id"
  )
  expect_named(snvs, expected_columns)
})
