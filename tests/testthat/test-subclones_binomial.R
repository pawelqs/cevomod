
test_that("predict_binoms works", {
  clones <- tibble(
    clone = c("Clone", "Subclone 1"),
    N_mutations = c(300, 100),
    cellularity = c(.5, .2)
  )
  binomial <- get_binomial_predictions(clones)
  map(binomial, sum)
  expect_equal(sum(binomial$Clone), 300)
  expect_equal(sum(binomial$`Subclone 1`), 100)
  expect_equal(max(binomial$Clone), 23.876771, tolerance = 0.0001)
  expect_equal(max(binomial$`Subclone 1`), 9.930021, tolerance = 0.0001)
})

