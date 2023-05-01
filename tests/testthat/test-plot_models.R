data("tcga_brca_test")

test_that("Comparing multiple models works", {
  object <- tcga_brca_test |>
    fit_powerlaw_tail_fixed(verbose = FALSE) |>
    fit_powerlaw_tail_optim(verbose = FALSE)
  model_names <- c("williams_neutral", "powerlaw_optim")
  column_name <- "powerlaw_pred"
  expect_warning({
    p <- compare_models(object, model_names, column_name)
  })
  expect_s3_class(p, c("gg", "ggplot"))
})


