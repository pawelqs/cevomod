set_cevomod_verbosity(0)

data("tcga_brca_test")
snvs <- SNVs(tcga_brca_test)
cnvs <- CNVs(tcga_brca_test)
meta <- tcga_brca_test$metadata |>
  mutate(
    patient_id = c("A", "A", "B", "B"),
    sample = c("S1", "S2", "S1", "S2")
  )


test_that("Filtering cevodata works", {
  patients <- unique(meta$patient_id)
  cd <- init_cevodata("TCGA BRCA small", cnvs = cnvs) |>
    add_CNV_data(cnvs = cnvs, "tcga2") |>
    add_SNV_data(snvs = snvs, "tcga") |>
    add_sample_data(meta) |>
    intervalize_mutation_frequencies() |>
    calc_SFS() |>
    calc_Mf_1f()
  cd$misc_by_sample[["slot1"]] <- meta$sample_id |>
    set_names(meta$sample_id) |>
    map(~c(1, 2, 3))
  cd$misc_by_patient[["slot2"]] <- patients |>
    set_names(patients) |>
    map(~c(1, 2, 3))
  cd <- filter(cd, sample_id == "TCGA-AC-A23H-01")
  expect_equal(nrow(cd$CNVs$cnvs), 285)
  expect_equal(nrow(cd$CNVs$tcga2), 285)
  expect_equal(nrow(cd$SNVs$tcga), 6419)
  expect_equal(unique(cd$CNVs$cnvs$sample_id), "TCGA-AC-A23H-01")
  expect_equal(unique(cd$CNVs$tcga2$sample_id), "TCGA-AC-A23H-01")
  expect_equal(unique(cd$SNVs$tcga$sample_id), "TCGA-AC-A23H-01")
  expect_equal(unique(cd$models$SFS$sample_id), "TCGA-AC-A23H-01")
  expected_misc_by_sample <- list(slot1 = list(`TCGA-AC-A23H-01` = c(1, 2, 3)))
  expect_equal(cd$misc_by_sample, expected_misc_by_sample)
  expected_misc_by_patient <- list(slot2 = list(A = c(1, 2, 3)))
  expect_equal(cd$misc_by_patient, expected_misc_by_patient)
})


test_that("Merging cevodata works", {
  cd <- init_cevodata("TCGA BRCA small", cnvs = cnvs) |>
    add_CNV_data(cnvs = cnvs, "tcga2") |>
    add_SNV_data(snvs = snvs, "tcga")
  cd1 <- filter(cd, sample_id == "TCGA-AC-A23H-01")
  cd2 <- filter(cd, sample_id != "TCGA-AC-A23H-01")
  cd_merged <- merge(cd1, cd2, name = "TCGA BRCA small", verbose = FALSE)
  output <- print(cd_merged) |>
    capture.output()
  expect_equal(object.size(cd_merged), object.size(cd))
})
