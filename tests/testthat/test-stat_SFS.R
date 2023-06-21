data("test_data")



test_that("Calculation of SFS works", {
  object <- test_data
  cd <- calc_SFS(object)
  expected <- read_tsv("../testdata/test_data_SFS.tsv", col_types = "ccdid")
  class(expected) <- c("cevo_SFS_tbl", class(expected))
  # write_tsv(cd$models$SFS, "tests/testdata/test_data_SFS.tsv")
  expect_equal(cd$models$SFS, expected)
})


test_that("plot(calc_SFS()) works", {
  dt <- SNVs(test_data)
  expect_warning({
    p <- dt %>%
      calc_SFS() %>%
      plot()
  })
  expect_s3_class(p, c("gg", "ggplot"))
  vdiffr::expect_doppelganger("test_data_sfs_1", p)
})


test_that("plot_SFS() works", {
  expect_warning({
    p <- test_data |>
      calc_SFS() |>
      plot_SFS()
  })
  vdiffr::expect_doppelganger("test_data_sfs_2", p)
})


test_that("get_non_zero_SFS_range works", {
  SFS <- tibble(
    sample_id = "S1",
    f = 1:100,
    y = 10
  ) |>
    mutate(
      y = case_when(
        f < 12 ~ 0,
        f == 40 ~ 0,
        f %in% c(54, 55) ~ 1,
        f > 65 ~ 0,
        TRUE ~ y
      ),
      f = f / 100
    )
  allowed_zero_bins <- 1
  y_treshold <- 1

  expected <- tibble(
    sample_id = "S1",
    from = 0.12,
    to = 0.65
  )
  res <- get_non_zero_SFS_range(SFS)

  expect_identical(res, expected)
})
