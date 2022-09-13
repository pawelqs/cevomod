test_that("%not in% works", {
  expect_true(1 %not in% c(0, 2, 3, 4))
  expect_false(1 %not in% c(0, 1, 2, 3, 4))
})


test_that("complete_missing_VAF_levels works", {
  dt <- tibble(
      sample_id = c("S1", "S1", "S2", "S2"),
      VAF = c(0.1, 0.2, 0.3, 0.4),
      y = 1
    ) |>
    group_by(sample_id)
  res <- complete_missing_VAF_levels(dt, fill = list(y = 0), digits = 1)
  expect_equal(min(res$VAF), 0)
  expect_equal(max(res$VAF), 1)
  expect_equal(nrow(res), 22)
})


test_that("complete_missing_VAF_levels works when df is groupped by two variables", {
  dt <- tibble(
      patient_id = c("S", "S", "S", "S", "A"),
      sample_id = c("S1", "S1", "S2", "S2", "A1"),
      VAF = c(0.1, 0.2, 0.3, 0.4, 0.4),
      y = 1
    ) |>
    group_by(patient_id, sample_id)
  res <- complete_missing_VAF_levels(dt, fill = list(y = 0), digits = 2)
  expect_equal(min(res$VAF), 0)
  expect_equal(max(res$VAF), 1)
  expect_equal(nrow(res), 303)
})


test_that("join_aes works", {
  aes_default <- aes(VAF, y, color = sample_id)
  aes_custom <- aes(color = patient_id, shape = effect)
  res <- join_aes(aes_default, aes_custom)
  expected <- aes(VAF, y, color = sample_id, shape = effect)
  expect_identical(res, expected)
})
