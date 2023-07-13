data("tcga_brca_test")
set_cevomod_verbosity(0)


test_that("Comparing multiple models works", {
  object <- tcga_brca_test |>
    fit_powerlaw_tail_fixed() |>
    fit_powerlaw_tail_optim()
  model_names <- c("powerlaw_fixed", "powerlaw_optim")
  column_name <- "powerlaw_pred"
  expect_warning({
    p <- compare_models(object, model_names, column_name)
  })
  expect_s3_class(p, c("gg", "ggplot"))
})


