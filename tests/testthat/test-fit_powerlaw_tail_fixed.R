data("tcga_brca_fitted")
verbose::verbose(cevoverse = 0)


# ---------------------------------- Fit ---------------------------------------

test_that("Fitting neutral partial models works", {
  rsq_treshold <- 0.98
  lm_length <- 0.05
  name <- "powerlaw_fixed"
  pct_left <- 0.05
  pct_right <- 0.95
  verbose <- get_verbosity()

  snvs <- SNVs(tcga_brca_fitted) |>
    filter(sample_id %in% c("TCGA-AC-A23H-01","TCGA-AN-A046-01"))
  object <- init_cevodata("Test") |>
    add_SNV_data(snvs) |>
    intervalize_mutation_frequencies() |>
    calc_Mf_1f() |>
    calc_SFS() |>
    fit_powerlaw_tail_fixed()
  models <- get_models(object)

  expected_coefs <- test_path("testdata", "tcga_brca_2samples.coefs_powerlaw_fixed.tsv") |>
    read_tsv(col_types = "cccdddddddl") |>
    mutate(
      model = "powerlaw_fixed",
      component = "Neutral tail",
      .before = "from"
    )
  expected_residuals <- test_path("testdata", "tcga_brca_2samples.residuals_powerlaw_fixed.tsv") |>
    read_tsv(show_col_types = FALSE)
  attr(expected_residuals, "f_column") <- "VAF"

  expect_s3_class(models, c("cv_powerlaw_models", "cv_subitem", "list"))
  expect_equal(models$coefs, expected_coefs)
  expect_equal(models$residuals, expected_residuals)
  expect_equal(models$info, list(f_column = "VAF"))
})


# ## Verify the above model visually to evaluate the changes
# powerlaw_fixed_old <- test_path("tcga_brca_partial_neutral_models.tsv") |>
#   read_tsv(col_types = "cccdddddddl") |>
#   mutate(
#     model = "powerlaw_fixed",
#     component = "Neutral tail",
#     .before = "from"
#   )
# object$models$powerlaw_fixed_old <- powerlaw_fixed_old
# object <- calc_powerlaw_model_residuals(object, models_name = "powerlaw_fixed_old")
# plot_models(object) +
#   geom_line(
#     aes(.data$f, .data$powerlaw_pred),
#     data = object$misc$residuals_powerlaw_fixed_old |>
#       filter(powerlaw_pred < 600),
#     color = "red", alpha = 0.5,
#     linetype = "dashed"
#   )



test_that("calc_powerlaw_curve() works", {
  dt <- tibble(VAF = 1:100/100, A = 176.5929, alpha = 2, nbins = 100)
  curve <- dt |>
    mutate(powerlaw_pred = calc_powerlaw_curve(VAF, A, alpha, nbins))
  expect_equal(nrow(curve), 100)
  expect_equal(
    curve$powerlaw_pred[1:5],
    c(17085.6315, 4271.4079, 1898.4035, 1067.8520, 683.4253),
    tolerance = 0.1
  )
  expect_equal(
    curve$powerlaw_pred[95:100],
    c(1.956708, 1.916156, 1.876851, 1.838743, 1.801784, 1.765929),
    tolerance = 0.0001
  )
  expect_true("powerlaw_pred" %in% names(curve))
})


object <- init_cevodata("Test") |>
  add_SNV_data(generate_neutral_snvs()) |>
  intervalize_mutation_frequencies() |>
  calc_Mf_1f() |>
  calc_SFS(bins = 100) |>
  fit_powerlaw_tail_fixed()


test_that("calc_powerlaw_model_residuals() creates proper tibble", {
  coefs <- get_models(object)$coefs |>
    filter(best)
  sfs <- get_SFS(object)
  resids <- calc_powerlaw_model_residuals(coefs, sfs)
  # resids <- get_residuals(cd, "powerlaw_fixed")
  expect_equal(nrow(resids), 101)
  expect_true(all(c("powerlaw_resid", "sampling_rate") %in% names(resids)))
  expect_equal(
    resids$sampling_rate[1:5],
    c(NaN, 1.000000, 0.999868, 0.99890, 0.995713),
    tolerance = 0.0001
  )
  expect_equal(is.na(resids) |> sum(), 1)
})


# ------------------------------- Plot -----------------------------------------

test_that("plot_Mf_1f_fits() works", {
  p <- plot_Mf_1f_fits(tcga_brca_fitted)
  vdiffr::expect_doppelganger("plot-Mf-1f-fits", p)
})


test_that("plot_neutral_A_coefficients() works", {
  p <- plot_neutral_A_coefficients(tcga_brca_fitted)
  vdiffr::expect_doppelganger("plot-neutral-A-coefficients", p)
})
