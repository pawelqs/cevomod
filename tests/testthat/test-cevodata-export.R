cd <- test_data
new_cnvs <- CNVs(cd) |>
  filter(sample_id != "Sample 3")
cd <- cd |>
  add_CNV_data(new_cnvs, "new_cnvs")
res <- to_clip(cd)


test_that("to_clip() creates correct tibbles for each sample", {
  expect_named(res, str_c("Sample ", 1:4))
  res |>
    map(\(x) names(x) == c("snvs", "cnvs", "purities")) |>
    map_lgl(all) |>
    all() |>
    expect_true()
})


test_that("to_clip() creates output files with correct column names", {
  expect_named(res$`Sample 1`$snvs, c("chromosome_index", "position", "alt_count", "ref_count"))
  expect_named(
    res$`Sample 1`$cnvs,
    c("chromosome_index", "start_position", "end_position", "major_cn", "minor_cn", "total_cn")
  )
  expect_type(res$`Sample 1`$purities, "double")
})
