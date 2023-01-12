data("tcga_brca_test")

test_that("Calculation of SFS works", {
  cd <- calc_SFS(tcga_brca_test)
  expected <- read_tsv("../testdata/tcga_brca_SFS.tsv", col_types = "ccdid")
  class(expected) <- c("cevo_SFS_tbl", class(expected))
  # write_tsv(cd$models$SFS, "tests/testdata/tcga_brca_SFS.tsv")
  expect_equal(cd$models$SFS, expected)
})


test_that("plot(calc_SFS()) works", {
  dt <- SNVs(tcga_brca_test)
  expect_warning({
    p <- dt %>%
      calc_SFS() %>%
      plot()
  })
  expect_s3_class(p, c("gg", "ggplot"))
  vdiffr::expect_doppelganger("plot(calc_SFS())", p)
})


test_that("plot_SFS() works", {
  expect_warning({
      p <- tcga_brca_test |>
    calc_SFS() |>
    plot_SFS()
  })
  vdiffr::expect_doppelganger("tcga_brca_test_SFS", p)
})
