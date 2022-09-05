data("tcga_brca_test")

test_that("Calculation of cumulative tails works", {
  dt <- tcga_brca_test %>%
    group_by(sample_id)
  res <- calc_cumulative_tails(dt)
  expected <- read_tsv("../testdata/tcga_brca_cumulative_tails.tsv", col_types = "cdiid")
  expect_identical(ungroup(res), expected)
})


test_that("Plotting of cumulative tails works", {
  dt <- tcga_brca_test %>%
    group_by(sample_id)
  p <- plot_cumulative_tails(dt)
  expect_s3_class(p, c("gg", "ggplot"))

  p2 <- dt %>%
    calc_cumulative_tails() %>%
    plot()
  expect_s3_class(p2, c("gg", "ggplot"))
})
