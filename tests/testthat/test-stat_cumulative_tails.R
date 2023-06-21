data("tcga_brca_test")

cd <- tcga_brca_test |>
  intervalize_mutation_frequencies(column = "VAF")

test_that("Calculation of cumulative tails works", {
  cd <- calc_cumulative_tails(cd)
  expected <- read_tsv("../testdata/tcga_brca_cumulative_tails.tsv", col_types = "ccdiid")
  # write_tsv(cd$models$cumulative_tails, "tests/testdata/tcga_brca_cumulative_tails.tsv")
  class(expected) <- c("cevo_cumulative_tails_tbl", class(expected))
  expect_equal(cd$models$cumulative_tails, expected)
})


test_that("Plotting of cumulative tails works", {
  p <- cd |>
    calc_cumulative_tails() |>
    plot_cumulative_tails()
  expect_s3_class(p, c("gg", "ggplot"))

  p2 <- SNVs(cd) %>%
    calc_cumulative_tails() %>%
    plot()
  expect_s3_class(p2, c("gg", "ggplot"))
})
