data("tcga_brca_test")

test_that("Calculation of Mf_1f works", {
  cd <- calc_Mf_1f(tcga_brca_test, verbose = FALSE)
  path <- test_path("tcga_brca_Mf_1f.tsv")
  expected <- read_tsv(path, col_types = "ccdiid")
  class(expected) <- c("cevo_Mf_1f_tbl", class(expected))
  # write_tsv(cd$models$Mf_1f, path)
  expect_equal(cd$models$Mf_1f, expected)
})


test_that("plot_Mf_1f() works", {
  p <- tcga_brca_test |>
    calc_Mf_1f(verbose = FALSE) |>
    plot_Mf_1f()
  vdiffr::expect_doppelganger("plot_Mf_1f", p)
})


test_that("plot(calc_Mf_1f()) works", {
  p <- SNVs(tcga_brca_test) |>
    calc_Mf_1f(verbose = FALSE) |>
    plot()
  vdiffr::expect_doppelganger("plot(calc_Mf_1f())", p)
})
