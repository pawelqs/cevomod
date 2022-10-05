data("tcga_brca_test")

test_that("Calculation of cumulative tails works", {
  cd <- calc_cumulative_tails(tcga_brca_test)
  expected <- read_tsv("../testdata/tcga_brca_cumulative_tails.tsv", col_types = "cdiid")
  # write_tsv(cd$models$cumulative_tails, "tests/testdata/tcga_brca_cumulative_tails.tsv")
  class(expected) <- c("cevo_cumulative_tails_tbl", class(expected))
  expect_identical(cd$models$cumulative_tails, expected)
})


test_that("Plotting of cumulative tails works", {
  p <- plot_cumulative_tails(tcga_brca_test)
  expect_s3_class(p, c("gg", "ggplot"))

  p2 <- SNVs(tcga_brca_test) %>%
    calc_cumulative_tails() %>%
    plot()
  expect_s3_class(p2, c("gg", "ggplot"))
})
