data("tcga_brca_test")

test_that("Fitting neutral partial models works", {
  snvs <- SNVs(tcga_brca_test) |>
    filter(sample_id %in% c("TCGA-AC-A23H-01","TCGA-AN-A046-01"))
  cd <- init_cevodata("Test") |>
    add_SNV_data(snvs) |>
    calc_Mf_1f() |>
    calc_SFS() |>
    fit_neutral_models(rsq_treshold = 0.99)
  expected <- "../testdata/tcga_brca_partial_neutral_models.tsv" |>
    read_tsv(col_types = "cdddddddl") |>
    mutate(
      model = "neutral_A/f^2",
      component = "Neutral tail",
      .before = "from"
    )
  class(expected) <- c("cevo_lm_models_tbl", class(expected))
  # write_tsv(cd$models$neutral_model, "tests/testdata/tcga_brca_partial_neutral_models.tsv")
  expect_equal(cd$models$neutral_model, expected)
})


cd <- init_cevodata("Test") |>
  add_SNV_data(generate_neutral_snvs()) |>
  calc_Mf_1f() |>
  calc_SFS() |>
  fit_neutral_models()


test_that("calc_powerlaw_curve works", {
  curve <- cd$models$neutral_model |>
    filter(best) |>
    calc_powerlaw_curve(binwidth = 0.01)
  expect_equal(nrow(curve), 100)
  expect_equal(
    curve$neutral_pred[1:5],
    c(17085.6315, 4271.4079, 1898.4035, 1067.8520, 683.4253),
    tolerance = 0.1
  )
  expect_equal(
    curve$neutral_pred[95:100],
    c(1.956708, 1.916156, 1.876851, 1.838743, 1.801784, 1.765929),
    tolerance = 0.0001
  )
  expect_true(all(c("neutr", "neutral_pred") %in% names(curve)))
})


test_that("calc_residuals creates proper tibble", {
  resids <- calc_residuals(cd)
  expect_equal(nrow(resids), 100)
  expect_true(all(c("neutral_resid", "sampling_rate") %in% names(resids)))
  expect_equal(
    resids$sampling_rate[1:5],
    c(1.000000, 0.999773, 0.998471, 0.994564, 0.987259),
    tolerance = 0.0001
  )
  expect_false(any(is.na(resids)))
})
