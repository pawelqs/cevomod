
# ------------------------- get_evolutionary_parameters ------------------------

test_that("get_evolutionary_parameters.cevodata() works", {
  Nmax = 10^10
  powerlaw_models <- tibble(
    sample_id = str_c("sample", 1:2),
    component = "powerlaw tail",
    A = c(100, 200),
    alpha = 2
  )
  binom <- tribble(
    ~sample_id, ~component,       ~N_mutations, ~frequency,
    "sample1",  "Clone",          100,          0.5,
    "sample1",  "Subclone 1",     100,          0.3,
    "sample2",  "Clone",          100,          0.5,
    "sample2",  "Subclone 1",     1000,         0.3,
    "sample2",  "Subclone 2",     1000,         0.1
  )
  models <- list(
    coefs = bind_rows(powerlaw_models, binom)
  )
  class(models) <- c("cv_powerlaw_subclones_models")
  object <- init_cevodata() |>
    add_models(models, "models")
  models_name <- active_models(object)

  res <- get_evolutionary_parameters(object)
  expect_equal(res$sample_id, c("sample1", "sample2", "sample2"))
  expect_equal(res$emergence_time, c(0.7213475, 3.6067376, 3.6067376), tolerance = 0.001)
  expect_equal(res$selection_coef, c(0.04190114, 0.19023756, 0.13216036), tolerance = 0.001)
})


test_that("get_evolutionary_parameters.cv_powerlaw_models() works", {
  models <- list(
    coefs = tibble(
      sample_id = str_c("sample", 0:5),
      component = "Neutral tail",
      A = 100,
      alpha = 2
    )
  )
  class(models) <- c("cv_powerlaw_models")

  res <- get_evolutionary_parameters(models)
  expect_equal(res$mutation_rate, rep(100, 6))
})


test_that("get_evolutionary_parameters.cv_powerlaw_subclones_models() works", {
  Nmax = 10^10
  powerlaw_models <- tibble(
    sample_id = str_c("sample", 1:2),
    component = "powerlaw tail",
    A = c(100, 200),
    alpha = 2
  )
  binom <- tribble(
    ~sample_id, ~component,       ~N_mutations, ~frequency,
    "sample1",  "Clone",          100,          0.5,
    "sample1",  "Subclone 1",     100,          0.3,
    "sample2",  "Clone",          100,          0.5,
    "sample2",  "Subclone 1",     1000,         0.3,
    "sample2",  "Subclone 2",     1000,         0.1
  )
  object <- list(
    coefs = bind_rows(powerlaw_models, binom)
  )
  class(object) <- c("cv_powerlaw_subclones_models")

  res <- get_evolutionary_parameters(object)
  expect_equal(res$sample_id, c("sample1", "sample2", "sample2"))
  expect_equal(res$emergence_time, c(0.7213475, 3.6067376, 3.6067376), tolerance = 0.001)
  expect_equal(res$selection_coef, c(0.04190114, 0.19023756, 0.13216036), tolerance = 0.001)
})


# ------------------------------- Other functions ------------------------------


test_that("mobster_emergence_time() works", {
  res <- mobster_emergence_time(100, 5)
  expect_equal(res, 14.4269504)
})


test_that("mobster_evolutionary_params() works", {
  subclones <- tribble(
    ~sample_id, ~component,    ~N_mutations, ~cellular_frequency, ~mutation_rate,
     "sample0",  "Subclone 1",  100,          0.1,                 100,
     "sample1",  "Subclone 1",  1000,         0.1,                 100,
     "sample2",  "Subclone 1",  1000,         0.1,                 200,
     "sample3",  "Subclone 1",  2000,         0.1,                 100,
     "sample4",  "Subclone 1",  2000,         0.3,                 100,
     "sample4",  "Subclone 2",  1000,         0.2,                 100,
     "sample5",  "Subclone 1",  2000,         0.7,                 100,
     "sample5",  "Subclone 2",  1000,         0.4,                 100,
  ) |>
    mutate(
      emergence_time = mobster_emergence_time(.data$N_mutations, .data$mutation_rate)
    ) |>
    nest(subclones = c("component", "N_mutations", "cellular_frequency", "emergence_time"))

  res <- subclones |>
    rowwise("sample_id", "mutation_rate") |>
    reframe(mobster_evolutionary_params(subclones))
  # write_tsv(res, "tests/testdata/expected_mobster_evolutionary_parameters.tsv")
  expected <- read_tsv(test_path("expected_mobster_evolutionary_parameters.tsv"), show_col_types = FALSE)
  expect_equal(res, expected)
})
