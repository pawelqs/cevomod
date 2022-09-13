test_that("Data generator generates correct data", {
  snvs <- generate_neutral_snvs()
  expect_s3_class(snvs, "tbl_df")
  expect_equal(nrow(snvs), 1532)
  expect_named(snvs, c("sample_id", "VAF", "n", "mut_id"))
})
