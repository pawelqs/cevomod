## code to prepare `snvs_tcga_brca` dataset goes here

# Go here -> https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
# Download the data using the Download Icon next to the dataset name
# Unpack
# Fix the path below

library(tidyverse)
set.seed(1234)


## ---------------------------------- SNVs ------------------------------------

mutations <- read_tsv("/mnt/dane/data/cbioportal/brca_tcga_pan_can_atlas_2018/data_mutations.txt")
glimpse(mutations)


snvs <- mutations %>%
  left_join(variant_classification) %>%
  transmute(
    sample_id = Tumor_Sample_Barcode,
    chrom = str_c("chr", Chromosome),
    pos = Start_Position,
    gene_symbol = Hugo_Symbol,
    ref = Reference_Allele,
    alt = Tumor_Seq_Allele2,
    ref_reads = t_ref_count,
    alt_reads = t_alt_count,
    VAF = alt_reads / (alt_reads + ref_reads),
    DP = ref_reads + alt_reads,
    ref_counts_normal = n_ref_count,
    alt_counts_normal = n_alt_count,
    variant_class = VARIANT_CLASS,
    variant_type = Variant_Type,
    Variant_Classification,
    variant_classification,
    impact = IMPACT,
    consequence = Consequence,
    NCBI_Build
  )
# class(snvs) <- c("cevo_snvs", class(snvs))

# usethis::use_data(snvs_tcga_brca, overwrite = TRUE)

## -------------------------------- CNVs --------------------------------------

cna_hg19 <- read_tsv("/mnt/dane/data/cbioportal/brca_tcga_pan_can_atlas_2018/data_cna_hg19.seg")
cnvs <- cna_hg19 |>
  transmute(
    sample_id = ID,
    chrom = str_c("chr", chrom),
    start = loc.start,
    end = loc.end,
    log_ratio = NA_real_,
    BAF = NA_real_,
    total_cn = NA_real_,
    major_cn = NA_real_,
    minor_cn = NA_real_,
    seg_mean = seg.mean
  )

## ---------------------------------- meta ------------------------------------

samples_data <- tibble(
  sample_id = unique(snvs$sample_id),
  patient_id = str_sub(sample_id, 1, 12),
  sample = "tumor"
)


clinical <- read_tsv(
  "/mnt/dane/data/cbioportal/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
  comment = "#"
) |>
  mutate(sex = str_to_lower(SEX)) |>
  select(patient_id = PATIENT_ID, sex) |>
  filter(patient_id %in% samples_data$patient_id)

TMB <- snvs %>%
  group_by(sample_id) %>%
  summarise(TMB = n()) %>%
  arrange(desc(TMB))

TMB |>
  mutate(
    sample_id = parse_factor(sample_id, levels = sample_id),
    sample_index = row_number()
  ) |>
  ggplot(aes(sample_index, TMB)) +
  geom_bar(stat = "identity")

## -------------------------------- cevodata ----------------------------------

tcga_brca <- init_cevodata("TCGA-BRCA data", genome = "hg37") |>
  add_SNV_data(snvs, name = "TCGA") |>
  add_CNV_data(cnvs, name = "TCGA") |>
  add_sample_data(samples_data) |>
  add_sample_data(TMB) |>
  add_patient_data(clinical) |>
  filter(TMB > 200)

top_mutated_patients <- TMB %>%
  pull(sample_id) %>%
  head(5) |>
  setdiff(c("TCGA-BH-A18G-01"))

tcga_brca_test <- tcga_brca |>
  filter(sample_id %in% top_mutated_patients) |>
  # calc_mutation_frequencies(method = "use_VAF") |>
  intervalize_mutation_frequencies() |>
  calc_SFS() |>
  calc_cumulative_tails() |>
  calc_Mf_1f() |>
  fit_powerlaw_tail_fixed() |>
  fit_subclones() |>
  fit_powerlaw_tail_optim() |>
  fit_subclones()
tcga_brca_test$active_models <- "powerlaw_fixed_subclones"


# usethis::use_data(tcga_brca, overwrite = TRUE)
usethis::use_data(tcga_brca_test, overwrite = TRUE)


