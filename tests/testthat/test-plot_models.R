data("tcga_brca_test")
set_cevomod_verbosity(0)


test_that("plot_models() works for powerlaw_optim_subclones", {
  models_name <- "powerlaw_optim_subclones"
  show_neutral_tail <- TRUE
  show_binomial_layer <- FALSE
  show_subclones <- TRUE
  show_final_fit <- TRUE
  nrow <- NULL
  ncol <- NULL

  suppressWarnings({
    p <- plot_models(test_data_fitted, models_name = "powerlaw_optim_subclones")
  })
  expect_s3_class(p, c("gg", "ggplot"))
  vdiffr::expect_doppelganger("powerlaw_optim_subclones", p)
})



test_that("plot_models() works for bs_powerlaw_model_bootstraps", {
  models_name <- "bs_powerlaw_model_bootstraps"
  show_neutral_tail <- TRUE
  show_binomial_layer <- FALSE
  show_subclones <- TRUE
  show_final_fit <- TRUE
  nrow <- NULL
  ncol <- NULL

  suppressWarnings({
    p <- plot_models(test_data_fitted, models_name = "bs_powerlaw_model_bootstraps")
  })
  expect_s3_class(p, c("gg", "ggplot"))
  vdiffr::expect_doppelganger("bs_powerlaw_model_bootstraps", p)
})



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
