
test_that("CliP export works", {
  cd <- test_data
  new_cnvs <- CNVs(cd) |>
    filter(sample_id != "Sample 3")
  cd <- cd |>
    add_CNV_data(new_cnvs, "new_cnvs")

  res <- to_clip(cd)

  expect_named(res, str_c("Sample ", 1:4))
  res |>
    map_lgl(\(x) all(names(x) == c("snvs", "cnvs", "purities"))) |>
    all() |>
    expect_true()
  res |>
    map("snvs") |>
    map_lgl(is_tibble) |>
    all() |>
    expect_true()
  res |>
    map("snvs") |>
    map_lgl(is_tibble) |>
    all() |>
    expect_true()
  res |>
    map("purities") |>
    map_lgl(is.double) |>
    all() |>
    expect_true()
})

