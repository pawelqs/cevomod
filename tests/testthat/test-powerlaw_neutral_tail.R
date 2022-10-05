cd <- init_cevodata("Test") |>
  add_SNV_data(generate_neutral_snvs()) |>
  calc_Mf_1f() |>
  calc_SFS() |>
  fit_neutral_lm()


test_that("calc_powerlaw_curve works", {
  curve <- cd$models$neutral_lm |>
    filter(best) |>
    calc_powerlaw_curve(binwidth = 0.01)
  expect_equal(nrow(curve), 100)
  expect_equal(
    curve$neutral_pred[1:5],
    c(17085.6315, 4271.4079, 1898.4035, 1067.8520, 683.4253),
    tolerance = 0.1
  )
  expect_equal(
    curve$neutral_pred[95:100],
    c(1.893145, 1.853910, 1.815882, 1.779012, 1.743254, 1.708563),
    tolerance = 0.0001
  )
  expect_true(all(c("neutr", "neutral_pred") %in% names(curve)))
})


test_that("layer_neutral_tail returns list of geoms", {
  geoms <- layer_neutral_tail(cd)
  expect_type(geoms, "list")
  expect_identical(map_chr(geoms, typeof) |> unique(), "environment")
})


test_that("calc_residuals creates proper tibble", {
  resids <- calc_residuals(cd)$models$residuals
  expect_equal(nrow(resids), 100)
  expect_true(all(c("neutral_resid", "sampling_rate") %in% names(resids)))
  expect_equal(
    resids$sampling_rate[1:5],
    c(1.0000000, 0.9997659, 0.9984197, 0.9943812, 0.9868310),
    tolerance = 0.0001
  )
  expect_false(any(is.na(resids)))
})
