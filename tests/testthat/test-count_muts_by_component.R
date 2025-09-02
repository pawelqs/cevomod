test_that("count_mutations_by_component() does not return Missing muts if nbins varies", {
  object <- test_data_fitted
  active_models(object) <- "powerlaw_optim_subclones"
  models_name <- active_models(object)

  expect_warning({
    res <- count_mutations_by_component(object, include_missing = TRUE)
  })
  res <- res |>
    filter(sample_id == "Sample 3")
  expected <- tibble(
    sample_id = "Sample 3",
    component = c("Clone", "Subclone 1", "Powerlaw tail"),
    N_mutations = c(712, 641, 2623),
    N = 3976
  ) |>
    mutate(frac = N_mutations / N)
  expect_equal(res, expected)
})


test_that("count_mutations_by_component() returns missing muts in nbins const.", {
  object <- test_data_fitted |>
    filter(sample_id == "Sample 3") |>
    calc_SFS(bins = 100) |>
    fit_powerlaw_tail_optim() |>
    fit_subclones(method = "mclust")
  models_name <- active_models(object)
  # plot_models(object)

  res <- count_mutations_by_component(object, include_missing = TRUE)
  expected <- tibble(
    sample_id = "Sample 3",
    component = c("Clone", "Subclone 1", "Powerlaw tail", "Missing mutations"),
    N_mutations = c(703, 924, 2320, 4931),
    N = 8878
  ) |>
    mutate(frac = N_mutations / N)
  expect_equal(res, expected)
})


test_that("count_powerlaw_tail_mutations() works", {
  object <- test_data_fitted
  models_name <- active_models(object)

  res <- count_powerlaw_tail_mutations(object)
  expected <- tibble(
    sample_id = str_c("Sample ", 1:4),
    component = "Powerlaw tail",
    N = c(1623, 3260, 2623, 2564)
  )
  expect_equal(res, expected)
})


test_that("count_missing_powerlaw_tail_mutations() works", {
  object <- test_data_fitted |>
    calc_SFS(bins = 100) |>
    fit_powerlaw_tail_optim()
  # plot_models(object)
  models_name <- active_models(object)
  min_f <- 0.01

  res <- count_missing_powerlaw_tail_mutations(object)
  expected <- tibble(
    sample_id = str_c("Sample ", 1:4),
    component = "Missing mutations",
    N = c(1179, 18695, 4931, 3534)
  )
  expect_equal(res, expected)
})



test_that("count_missing_powerlaw_tail_mutations() throws warning if bins number varies", {
  object <- test_data_fitted
  active_models(object) <- "powerlaw_optim_subclones"
  models_name <- active_models(object)
  min_f <- 0.01

  expect_warning({
    res <- count_missing_powerlaw_tail_mutations(object)
  })
  expect_null(res)
})

