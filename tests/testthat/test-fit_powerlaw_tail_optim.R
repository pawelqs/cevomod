data("tcga_brca_fitted")
verbose::verbose(cevoverse = 0)


# --------------------------- fit_powerlaw_tail_optim -------------------------

test_that("fit_powerlaw_tail_optim() returns models with non-negative objective values", {
  name <- "powerlaw_optim"
  bootstraps <- FALSE
  allowed_zero_bins <- 2
  y_treshold <- 1
  y_threshold_pct <- 0.01
  av_filter <- c(1/3, 1/3, 1/3)
  peak_detection_upper_limit <- 0.3
  reward_upper_limit <- 0.4
  control <- list(maxit = 1000, ndeps = c(0.1, 0.01))
  verbose <- get_verbosity()

  object <- test_data_fitted |>
    intervalize_mutation_frequencies() |>
    fit_powerlaw_tail_optim()
  td <- get_model_coefficients(object, "powerlaw_optim")

  expected_coefs <- test_path("testdata", "test_data.coefs_powerlaw_optim.tsv") |>
    read_tsv(show_col_types = FALSE)
  expect_equal(td, expected_coefs)


  expected <- test_path("testdata", "test_data.residuals_powerlaw_optim.tsv") |>
    read_tsv(show_col_types = FALSE)
  attr(expected, "f_column") <- "VAF"
  expect_equal(get_model_residuals(object), expected)

  expect_s3_class(get_models(object, "powerlaw_optim"), "cv_powerlaw_models")
})


test_that("fit_powerlaw_tail_optim() bootstrapping test", {
  name <- "powerlaw_optim"
  bootstraps <- 2
  allowed_zero_bins <- 2
  y_treshold <- 1
  y_threshold_pct <- 0.01
  av_filter <- c(1/3, 1/3, 1/3)
  peak_detection_upper_limit <- 0.3
  reward_upper_limit <- 0.4
  control <- list(maxit = 1000, ndeps = c(0.1, 0.01))
  verbose <- get_verbosity()
  object <- tcga_brca_fitted |>
    intervalize_mutation_frequencies() |>
    filter(sample_id %in% tcga_brca_fitted$metadata$sample_id[1:2])

  suppressWarnings({
    object <- fit_powerlaw_tail_optim(object, bootstraps = bootstraps)
  })

  models <- get_models(object)
  expect_s3_class(models, "cv_powerlaw_models")

  # Models
  expected_columns <- c(
    "sample_id", "model", "component",
    "A", "A.lower", "A.upper", "alpha", "alpha.lower", "alpha.upper"
  )
  expect_named(models$coefs, expected_columns)
  expect_equal(nrow(models$coefs), 2)
  # Bootstrap coefs
  expect_equal(nrow(models$bootstrap_coefs), 4)
  expect_true(all(models$bootstrap_models$value > 1990))
  # Residuals
  expect_equal(nrow(models$residuals), 176)
  # Info
  expect_equal(models$info, list(f_column = "VAF"))
})


# ---------------------- fit_powerlaw_tail_optim methods -----------------------


test_that("fit_powerlaw_tail_optim.cevo_SFS_tbl() works", {
  name <- "powerlaw_optim"
  bootstraps <- FALSE
  allowed_zero_bins <- 2
  y_treshold <- 1
  y_threshold_pct <- 0.01
  av_filter <- c(1/3, 1/3, 1/3)
  peak_detection_upper_limit <- 0.3
  reward_upper_limit <- 0.4
  control <- list(maxit = 1000, ndeps = c(0.1, 0.01))
  verbose <- get_verbosity()

  object <- tcga_brca_fitted |>
    intervalize_mutation_frequencies() |>
    get_SFS()

  res <- fit_powerlaw_tail_optim(object)
  expect_s3_class(res, c("cv_powerlaw_models", "list"))
  # expect_equal(res$coefs, 2)
  expect_true(all(res$coefs$value > 1990))
})


test_that("fit_powerlaw_tail_optim.cevo_SFS_bootstraps() works", {
  name <- "powerlaw_optim"
  bootstraps <- 2
  allowed_zero_bins <- 2
  y_treshold <- 1
  y_threshold_pct <- 0.01
  av_filter <- c(1/3, 1/3, 1/3)
  peak_detection_upper_limit <- 0.3
  reward_upper_limit <- 0.4
  control <- list(maxit = 1000, ndeps = c(0.1, 0.01))
  verbose <- TRUE

  object <- tcga_brca_fitted |>
    intervalize_mutation_frequencies() |>
    filter(sample_id %in% tcga_brca_fitted$metadata$sample_id[1:2]) |>
    calc_SFS_resamples(times = bootstraps)
  object <- object$`TCGA-AC-A23H-01`

  suppressWarnings({
    res <- fit_powerlaw_tail_optim(object, bootstraps = bootstraps)
  })

  expect_s3_class(res, "cv_powerlaw_models")

  expected_columns <- c(
    "sample_id", "model", "component",
    "A", "A.lower", "A.upper", "alpha", "alpha.lower", "alpha.upper"
  )
  expect_named(res$coefs, expected_columns)
  expect_equal(nrow(res$coefs), 1)

  expect_equal(nrow(res$bootstrap_coefs), 2)
  expect_true(all(res$bootstrap_coefs$value > 1990))
})


# --------------------------------- Helpers ------------------------------------

test_that("td_objective_function() works", {
  object <- tcga_brca_fitted
  sfs <- get_SFS(object)
  bounds <- get_f_range(SNVs(object), pct_left = 0.02, pct_right = 0.98)
  nbins <- SNVs(object) |>
    group_by(.data$sample_id) |>
    summarise(nbins = median(DP))

  data <- sfs |>
    left_join(bounds, by = "sample_id") |>
    filter(.data$f > .data$lower_bound, .data$f < .data$higher_bound) |>
    select("sample_id", "f", "y") |>
    nest_by(.data$sample_id) |>
    left_join(nbins, by = "sample_id")

  dt <- data$data[[1]]
  x <- dt$f
  y <- dt$y
  params <- c(255.9222929 / 74, 2.912159)
  res <- -td_objective_function(params, x, y)
  expect_true(res > 2000)
})


test_that("get_non_zero_SFS_range() works", {
  SFS <- tibble(
    sample_id = "S1",
    f = 1:100,
    y = 10
  ) |>
    mutate(
      y = case_when(
        f < 12 ~ 0,
        f == 40 ~ 0,
        f %in% c(54, 55) ~ 1,
        f > 65 ~ 0,
        TRUE ~ y
      ),
      f = f / 100
    )
  allowed_zero_bins <- 1
  y_treshold <- 1

  expected <- tibble(
    sample_id = "S1",
    from = 0.12,
    to = 0.65
  )
  res <- get_non_zero_SFS_range(SFS)

  expect_identical(res, expected)
})


test_that("Calculation of SFS with resampling works", {
  cd <- cevodata::test_data
  sfs_resamples <- calc_SFS_resamples(cd, times = 2)

  expect_type(sfs_resamples, "list")
  expect_equal(length(sfs_resamples), 4)

  sfs_resamples |>
    walk(expect_s3_class, "cevo_SFS_bootstraps")

  attribs <- attributes(sfs_resamples$`Sample 1`$sfs[[1]])
  expect_equal(attribs$f_column, "VAF")
})


test_that("Calculation of SFS with resampling works with CCF/2", {
  cd <- cevodata::test_data |>
    calc_mutation_frequencies()
  sfs_resamples <- calc_SFS_resamples(cd, times = 2)

  expect_type(sfs_resamples, "list")
  expect_equal(length(sfs_resamples), 4)

  sfs_resamples |>
    walk(expect_s3_class, "cevo_SFS_bootstraps")

  attribs <- attributes(sfs_resamples$`Sample 1`$sfs[[1]])
  expect_equal(attribs$f_column, "CCF/2")
})

