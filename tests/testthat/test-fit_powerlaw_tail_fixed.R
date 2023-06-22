data("tcga_brca_test")
set_cevomod_verbosity(0)


test_that("Fitting neutral partial models works", {
  snvs <- SNVs(tcga_brca_test) |>
    filter(sample_id %in% c("TCGA-AC-A23H-01","TCGA-AN-A046-01"))
  object <- init_cevodata("Test") |>
    add_SNV_data(snvs) |>
    intervalize_mutation_frequencies() |>
    calc_Mf_1f() |>
    calc_SFS() |>
    fit_powerlaw_tail_fixed(rsq_treshold = 0.99)
  expected <- "../testdata/tcga_brca_partial_neutral_models.tsv" |>
    read_tsv(col_types = "cccdddddddl") |>
    mutate(
      model = "powerlaw_fixed",
      component = "Neutral tail",
      .before = "from"
    )
  class(expected) <- c("cevo_powerlaw_models", class(expected))
  # write_tsv(object$models$powerlaw_fixed, "tests/testdata/tcga_brca_partial_neutral_models.tsv")
  expect_equal(get_models(object, best_only = FALSE), expected)
})


# ## Verify the above model visually to evaluate the changes
# powerlaw_fixed_old <- "tests/testdata/tcga_brca_partial_neutral_models.tsv" |>
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
  cd <- calc_powerlaw_model_residuals(object, "powerlaw_fixed")
  resids <- get_residuals(cd, "powerlaw_fixed")
  expect_equal(nrow(resids), 101)
  expect_true(all(c("powerlaw_resid", "sampling_rate") %in% names(resids)))
  expect_equal(
    resids$sampling_rate[1:5],
    c(NaN, 1.000000, 0.999868, 0.99890, 0.995713),
    tolerance = 0.0001
  )
  expect_equal(is.na(resids) |> sum(), 1)
})

