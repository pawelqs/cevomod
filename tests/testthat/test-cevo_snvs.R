
test_that("as_cevo_snvs() works", {
  snvs <- tibble(sample_id = "A", chrom = "chr", pos = 1:3, VAF = 0.5)
  expected <- tibble(
    sample_id = "A",
    chrom = "chr",
    pos = 1:3,
    gene_symbol = NA_character_,
    ref = NA_character_,
    alt = NA_character_,
    ref_reads = NA_real_,
    alt_reads = NA_real_,
    impact = NA_character_,
    VAF = 0.5
  )
  class(expected) <- c("cevo_snvs", "tbl_df", "tbl", "data.frame")
  expect_identical(as_cevo_snvs(snvs), expected)
})
