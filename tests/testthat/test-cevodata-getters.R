
test_that("SNVs_CNVs() works", {
  snvs <- tribble(
    ~sample_id, ~chrom, ~pos, ~VAF,
    "S1",       "chr1", 100,  0.1,
    "S1",       "chr1", 200,  0.1,
    "S1",       "chr1", 300,  0.1,
    "S1",       "chr1", 400,  0.1,
    "S1",       "chr2", 100,  0.1,
    "S2",       "chr1", 100,  0.1,
    "S2",       "chr1", 300,  0.1
  )
  cnvs <- tribble(
    ~sample_id, ~chrom, ~start, ~end, ~total_cn, ~minor_cn, ~normal_cn,
    "S1",       "chr1", 1,      250,  3,         1,         2,
    "S1",       "chr1", 350,    450,  4,         1,         2,
    "S1",       "chr2", 1,      250,  2,         1,         2,
    "S2",       "chr1", 1,      250,  3,         1,         2
  )
  cd <- init_cevodata("test", snvs = snvs, cnvs = cnvs)
  expected <- tribble(
    ~sample_id, ~chrom, ~pos, ~VAF, ~total_cn, ~minor_cn, ~normal_cn,
    "S1",       "chr1", 100,  0.1,  3,         1,         2,
    "S1",       "chr1", 200,  0.1,  3,         1,         2,
    "S1",       "chr1", 300,  0.1,  NA_real_,  NA_real_,  NA_real_,
    "S1",       "chr1", 400,  0.1,  4,         1,         2,
    "S1",       "chr2", 100,  0.1,  2,         1,         2,
    "S2",       "chr1", 100,  0.1,  3,         1,         2,
    "S2",       "chr1", 300,  0.1,  NA_real_,  NA_real_,  NA_real_
  )
  res <- SNVs_CNVs(cd)
  # sample_data <- tibble(sample_id = c("S1", "S2"), sex = c())
  expect_identical(res$total_cn, expected$total_cn)
  expect_identical(res$minor_cn, expected$minor_cn)
})
