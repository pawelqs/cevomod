data("tcga_brca_test")

test_that("Calculation of SFS works", {
  dt <- tcga_brca_test %>%
    group_by(sample_id)
  res <- calc_SFS(dt)
  expected <- read_tsv("../testdata/tcga_brca_SFS.tsv", col_types = "cdid")
  expect_identical(ungroup(res), expected)
})


test_that("plot(calc_SFS()) works", {
  dt <- tcga_brca_test %>%
    group_by(sample_id)
  p <- dt %>%
    calc_SFS() %>%
    plot()
  expect_s3_class(p, c("gg", "ggplot"))
  vdiffr::expect_doppelganger("plot(calc_SFS())", p)
})


test_that("plot_SFS() works", {
  dt <- tcga_brca_test %>%
    group_by(sample_id)
  p <- plot_SFS(dt)
  vdiffr::expect_doppelganger("plot_SFS()", p)
})
