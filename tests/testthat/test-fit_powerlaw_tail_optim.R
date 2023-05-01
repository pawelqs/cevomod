data("tcga_brca_test")


test_that("fit_powerlaw_tail_optim models have non-negative objective fun value", {
  cd <- tcga_brca_test |>
    fit_powerlaw_tail_optim(verbose = FALSE)
  td <- get_models(cd, "powerlaw_optim")
  expect_true(all(td$value > 1990))
  expect_s3_class(get_powerlaw_models(cd, "powerlaw_optim"), "cevo_powerlaw_models")
})


test_that("Testing td_objective_function", {
  object <- tcga_brca_test
  sfs <- get_SFS(object)
  bounds <- get_VAF_range(SNVs(object), pct_left = 0.02, pct_right = 0.98)
  nbins <- get_sample_sequencing_depths(SNVs(object)) |>
    transmute(.data$sample_id, nbins = .data$median_DP)

  data <- sfs |>
    left_join(bounds, by = "sample_id") |>
    filter(.data$VAF > .data$lower_bound, .data$VAF < .data$higher_bound) |>
    select("sample_id", "VAF", "y") |>
    nest_by(.data$sample_id) |>
    left_join(nbins, by = "sample_id")

  dt <- data$data[[1]]
  x <- dt$VAF
  y <- dt$y
  params <- c(255.9222929 / 74, 2.912159)
  res <- -td_objective_function(params, x, y)
  expect_true(res > 2000)
})


# test_that("Testing tung-durret models on Shlush_AML and tcga_brca_test data", {
  # object <- tcga_brca_test |>
  #   fit_williams_neutral_models() |>
  #   fit_powerlaw_tail_optim()
  # model_names <- c("williams_neutral", "powerlaw_optim")
  # column_name <- "powerlaw_pred"
  # compare_models(object, model_names, column_name)
  #
  # data("Shlush_AML", package = "cevoDatasets")
  # cd <- Shlush_AML |>
  #   calc_SFS() |>
  #   calc_Mf_1f() |>
  #   fit_williams_neutral_models() |>
  #   fit_powerlaw_tail_optim()
  # compare_models(cd, model_names, column_name)
  # cd |>  filter(sample_id == "AMLRO-9_Rx") |> compare_models(model_names, column_name)
# })
