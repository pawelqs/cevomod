data("tcga_brca_test")

test_that("cut_f() properly names columns", {
  breaks <- SNVs(tcga_brca_test) |>
    get_interval_breaks()
  breaks <- breaks$`TCGA-AC-A23H-01`
  snvs <- SNVs(tcga_brca_test) |>
    filter(sample_id =="TCGA-AC-A23H-01")
  res <- cut_f(snvs, breaks, column = "VAF")
  expect_true("VAF_interval" %in% names(res))

  snvs <- tcga_brca_test |>
    calc_mutation_frequencies() |>
    SNVs() |>
    filter(sample_id =="TCGA-AC-A23H-01")
  res <- cut_f(snvs, breaks)
  expect_true("f_interval" %in% names(res))
})
