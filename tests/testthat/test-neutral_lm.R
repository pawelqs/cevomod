data("snvs_test")

test_that("Fitting neutral partial models works", {
  models <- snvs_test |>
    filter(sample_id %in% c("TCGA-AC-A23H-01","TCGA-AN-A046-01")) |>
    group_by(sample_id) |>
    calc_Mf_1f() |>
    fit_neutral_lm(rsq_treshold = 0.99)
  expected <- read_tsv("../testdata/tcga_brca_partial_neutral_models.tsv", col_types = "cddddddl")
  expect_equal(ungroup(models), expected)
})
