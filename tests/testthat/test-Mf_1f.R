data("tcga_brca_test")

test_that("Calculation of Mf_1f works", {
  dt <- tcga_brca_test %>%
    group_by(sample_id)
  res <- calc_Mf_1f(dt)
  expected <- read_tsv("../testdata/tcga_brca_Mf_1f.tsv", col_types = "cdiid")
  expect_identical(ungroup(res), expected)
})
