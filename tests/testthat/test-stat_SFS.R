data("tcga_brca_test")

test_that("Calculation of SFS works", {
  cd <- calc_SFS(tcga_brca_test)
  expected <- read_tsv("../testdata/tcga_brca_SFS.tsv", col_types = "cdid")
  # write_tsv(cd$models$SFS, "tests/testdata/tcga_brca_SFS.tsv")
  expect_identical(cd$models$SFS, expected)
})


test_that("plot(calc_SFS()) works", {
  dt <- SNVs(tcga_brca_test) %>%
    group_by(sample_id)
  p <- dt %>%
    calc_SFS() %>%
    plot()
  expect_s3_class(p, c("gg", "ggplot"))
  vdiffr::expect_doppelganger("plot(calc_SFS())", p)
})


test_that("plot_SFS() works", {
  p <- plot_SFS(tcga_brca_test)
  vdiffr::expect_doppelganger("plot_SFS()", p)
})
