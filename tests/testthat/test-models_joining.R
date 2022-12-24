data("AMLRO", package = "cevoDatasets")
object <- AMLRO |>
  filter(.data$patient_id %in% c("AMLRO-15")) |>
  run_cevomod()
# snvs <- SNVs(object)


test_that("get_interval_centers() works", {
  intervals <- (1:1000/1000) |>
    cut_width(width = 0.1, center = 0.1) |>
    levels()
  res <- get_interval_centers(intervals)
  expect_equal(res, 0:10/10)
})


test_that("get_interval_width() works", {
  intervals <- (1:1000/1000) |>
    cut_width(width = 0.1, center = 0.1) |>
    levels()
  res <- get_interval_width(intervals)
  expect_equal(res, 0.1)
})

