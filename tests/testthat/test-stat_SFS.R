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


test_that("get_non_zero_SFS_range works", {
  SFS <- tibble(
    sample_id = "S1",
    VAF = 1:100,
    y = 10
  ) |>
    mutate(
      y = case_when(
        VAF < 12 ~ 0,
        VAF == 40 ~ 0,
        VAF %in% c(54, 55) ~ 1,
        VAF > 65 ~ 0,
        TRUE ~ y
      ),
      VAF = VAF / 100
    )
  allowed_zero_bins <- 1
  y_treshold <- 1
  expected <- tibble(
    sample_id = "S1",
    from = 0.12,
    to = 0.65
  )
  expect_identical(get_non_zero_SFS_range(SFS), expected)
})
