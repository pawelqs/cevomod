test_that("mobster_evolutionary_parameters() works", {
  subclones <- tribble(
    ~sample_id, ~component,    ~N_mutations, ~subclone_frequency, ~mutation_rate_williams,
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
      emergence_time = get_emergence_time(.data$N_mutations, .data$mutation_rate_williams)
    ) |>
    nest(subclones = c("component", "N_mutations", "subclone_frequency", "emergence_time"))

  res <- subclones |>
    rowwise("sample_id") |>
    reframe(mobster_evolutionary_parameters(subclones, mu))
  # write_tsv(res, "tests/testdata/expected_mobster_evolutionary_parameters.tsv")
  expected <- read_tsv("../testdata/expected_mobster_evolutionary_parameters.tsv", show_col_types = FALSE)
  expect_equal(res, expected)
})
