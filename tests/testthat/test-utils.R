test_that("%not in% works", {
  expect_true(1 %not in% c(0, 2, 3, 4))
  expect_false(1 %not in% c(0, 1, 2, 3, 4))
})


test_that("join_aes works", {
  aes_default <- aes(VAF, y, color = sample_id)
  aes_custom <- aes(color = patient_id, shape = effect)
  res <- join_aes(aes_default, aes_custom)
  expected <- aes(VAF, y, color = patient_id, shape = effect)
  expect_identical(res, expected)
})
