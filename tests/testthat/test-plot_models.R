data("tcga_brca_test")

test_that("Comparing multiple models works", {
  object <- tcga_brca_test |>
    fit_neutral_models(verbose = FALSE) |>
    fit_tung_durrett_models(verbose = FALSE)
  model_names <- c("williams_neutral", "tung_durrett")
  column_name <- "powerlaw_pred"
  expect_warning({
    p <- compare_models(object, model_names, column_name)
  })
  expect_s3_class(p, c("gg", "ggplot"))
})


