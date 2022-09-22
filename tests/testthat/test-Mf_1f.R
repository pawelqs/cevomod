data("tcga_brca_test")

test_that("Calculation of Mf_1f works", {
  cd <- calc_Mf_1f(tcga_brca_test)
  expected <- read_tsv("../testdata/tcga_brca_Mf_1f.tsv", col_types = "cdiid")
  cd$models$Mf_1f |>
    select(-patient_id, -sample) |>
    expect_identical(expected)
  # expect_identical(cd$models$Mf_1f, expected)
})


test_that("plot_Mf_1f() works", {
  p <- tcga_brca_test |>
    plot_Mf_1f()
  vdiffr::expect_doppelganger("plot_Mf_1f", p)
})


test_that("plot(calc_Mf_1f()) works", {
  p <- SNVs(tcga_brca_test) |>
    group_by(sample_id) |>
    calc_Mf_1f() |>
    plot()
  vdiffr::expect_doppelganger("plot(calc_Mf_1f())", p)
})
