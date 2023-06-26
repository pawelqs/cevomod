set_cevomod_verbosity(0)


test_that("calc_mutation_frequencies() adds correct CCF column", {
  snvs <- tribble(
    ~sample_id, ~chrom, ~pos, ~VAF,
    "S1",       "chr1", 100,  0.1,
    "S1",       "chr1", 200,  0.1,
    "S1",       "chr1", 300,  0.9,  # no CNV, will be dropped
    "S1",       "chr1", 400,  0.4,
    "S1",       "chr2", 100,  0.2,
    "S2",       "chr1", 100,  0.16,
    "S2",       "chr1", 300,  0.16  # no CNV, will be dropped
  )
  cnvs <- tribble(
    ~sample_id, ~chrom, ~start, ~end, ~total_cn, ~minor_cn, ~normal_cn,
    "S1",       "chr1", 1,      250,  3,         1,         2,
    "S1",       "chr1", 350,    450,  3,         1,         2,
    "S1",       "chr2", 1,      250,  3,         1,         2,
    "S2",       "chr1", 1,      250,  3,         1,         2
  )
  purities <- tibble(
    sample_id = c("S1", "S2"),
    purity = c(0.5, 1)
  )
  object <- init_cevodata("test", snvs = snvs, cnvs = cnvs) |>
    add_sample_data(purities)

  expected <- object |>
    SNVs() |>
    mutate(
      CCF = c(0.5, 0.5, NA_real_, 1, 1, 0.48, NA_real_),
      `CCF/2` = CCF / 2
    )
  res <- object |>
    calc_mutation_frequencies(method = "Dentro", rm_intermediate_cols = TRUE)

  expect_identical(names(res$SNVs), names(object$SNVs))
  expect_equal(SNVs(res), expected)
})


test_that("dentro_2015_correction() works", {
  snvs <- tibble(
    VAF = c(0.1, 0.4, 0.2),
    total_cn = c(3, 3, 3),
    normal_cn = c(2, 2, 2),
    purity = 0.5
  )
  res <- dentro_2015_correction(snvs)
  expect_named(res, c("VAF", "CCF", "CCF/2", "total_cn", "normal_cn", "purity", "u", "m"))
  expect_equal(res$CCF, c(0.5, 1, 1))
})
