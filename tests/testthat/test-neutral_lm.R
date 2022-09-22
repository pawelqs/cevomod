data("tcga_brca_test")

test_that("Fitting neutral partial models works", {
  snvs <- SNVs(tcga_brca_test) |>
    filter(sample_id %in% c("TCGA-AC-A23H-01","TCGA-AN-A046-01"))
  cd <- init_cevodata("Test") |>
    add_SNV_data(snvs) |>
    calc_Mf_1f() |>
    fit_neutral_lm(rsq_treshold = 0.99)
  expected <- read_tsv("../testdata/tcga_brca_partial_neutral_models.tsv", col_types = "cccddddddl")
  class(expected) <- c("cevo_lm_models_tbl", class(expected))
  # write_tsv(cd$models$neutral_lm, "tests/testdata/tcga_brca_partial_neutral_models.tsv")
  expect_equal(cd$models$neutral_lm, expected)
})
