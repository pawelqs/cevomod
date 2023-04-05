test_that("filter_SNVs_by_regions works", {
  snvs <- tibble(
    sample_id = "S1",
    chrom = str_c("chr", c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)),
    pos = c(100:104, 100:104)
  )
  regions <- tibble(
    chrom = c("chr1", "chr2"),
    start = c(103, 100),
    end = c(200, 102)
  )
  bed_file = "../testdata/regions.tsv"
  expected <- tibble(
    sample_id = "S1",
    chrom = str_c("chr", c(1, 1, 2, 2, 2)),
    pos = c(103L, 104L, 100L, 101L, 102L)
  )

  expect_identical(filter_SNVs_by_regions(snvs, regions = regions), expected)
  expect_identical(filter_SNVs_by_regions(snvs, bed_file = bed_file), expected)
})
