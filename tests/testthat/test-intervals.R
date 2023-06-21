data("test_data")



test_that("cut_f_intervals() works correctly", {
  snvs <- test_data |>
    filter(sample_id =="Sample 1") |>
    SNVs()
  column <- "VAF"
  res <- cut_f_intervals(snvs, column = column, bins = 10)
  expect_true( all(c("f_interval", "f") %in% names(res)) )

  interval_centers <- as.character(unique(res$f))
  allowed_interval_centers <- as.character(seq(-5, 95, by = 5) / 100)
  expect_equal( length(setdiff(interval_centers, allowed_interval_centers)), 0 )
})

