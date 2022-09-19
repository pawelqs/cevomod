data("snvs_test")

test_that("Calculation of Mf_1f works", {
  dt <- snvs_test %>%
    group_by(sample_id)
  res <- calc_Mf_1f(dt)
  expected <- read_tsv("../testdata/tcga_brca_Mf_1f.tsv", col_types = "cdiid")
  expect_identical(ungroup(res), expected)
})


test_that("plot_Mf_1f() works", {
  p <- snvs_test |>
    group_by(sample_id) |>
    plot_Mf_1f()
  vdiffr::expect_doppelganger("plot_Mf_1f", p)
})


test_that("plot(calc_Mf_1f()) works", {
  p <- snvs_test |>
    group_by(sample_id) |>
    calc_Mf_1f() |>
    plot()
  vdiffr::expect_doppelganger("plot(calc_Mf_1f())", p)
})
