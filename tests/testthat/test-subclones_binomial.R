
test_that("predict_binoms works", {
  x <- predict_binoms(Ns = c(300, 100), means = c(.2, .5))$binom_pred
  # x <- predict_binoms(Ns = c(300, 100), means = c(.2, .5))
  expect_equal(sum(x), 400)
  expect_equal(max(x), 29.7900645)
})

