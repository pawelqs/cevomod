snvs <- generate_neutral_snvs() |>
  group_by(sample_id)
mf1f <- calc_Mf_1f(snvs)
sfs <- calc_SFS(snvs)
lm <- fit_neutral_lm(mf1f)

test_that("calc_powerlaw_curve works", {
  curve <- lm |>
    slice(1) |>
    calc_powerlaw_curve()
  expect_equal(nrow(curve), 100)
  expect_equal(curve$n[1:5], c(18984, 4746, 2109, 1187, 759), tolerance = 0.1)
  expect_equal(
    curve$n[95:100],
    c(2.103494, 2.059900, 2.017646, 1.976680, 1.936949, 1.898404),
    tolerance = 0.0001
  )
  expect_true(all(c("neutr", "n") %in% names(curve)))
})


test_that("layer_neutral_tail returns list of geoms", {
  geoms <- layer_neutral_tail(lm)
  expect_type(geoms, "list")
  expect_identical(map_chr(geoms, typeof) |> unique(), "environment")
})


test_that("estimate_sampling_rate returns proper tibble", {
  sampling_rate <- estimate_sampling_rate(sfs, lm)
  expect_equal(nrow(sampling_rate), 57)
  expect_true(all(c("err", "sampling_rate") %in% names(sampling_rate)))
  expect_equal(
    sampling_rate$sampling_rate[1:5],
    c(1.0000000, 0.9997893, 0.9985778, 0.9949431, 0.9881479),
    tolerance = 0.0001
  )
})
