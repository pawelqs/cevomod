test_that("Data generator generates correct data", {
  snvs <- generate_neutral_snvs()
  expect_s3_class(snvs, "tbl_df")
  expect_equal(nrow(snvs), 1532)
  expected_columns <- c(
    "patient_id", "sample_id", "VAF", "DP", "n", "alt_reads", "ref_reads", "mutation_id"
  )
  expect_named(snvs, expected_columns)
})
