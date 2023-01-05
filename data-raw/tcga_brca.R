## code to prepare `snvs_tcga_brca` dataset goes here

# Go here -> https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
# Download the data using the Download Icon next to the dataset name
# Unpack
# Fix the path below

library(tidyverse)

mutations <- read_tsv("/mnt/dane/data/brca_tcga_pan_can_atlas_2018/data_mutations.txt")
glimpse(mutations)


snvs_tcga_brca <- mutations %>%
  left_join(variant_classification) %>%
  transmute(
    sample_id = Tumor_Sample_Barcode,
    chrom = Chromosome,
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
    variant_type =  Variant_Type,
    Variant_Classification,
    variant_classification,
    impact = IMPACT,
    consequence = Consequence,
    NCBI_Build
  )
class(snvs_tcga_brca) <- c("cevo_snvs", class(snvs_tcga_brca))

# usethis::use_data(snvs_tcga_brca, overwrite = TRUE)

top_mutated_patients <- snvs_tcga_brca %>%
  group_by(sample_id) %>%
  count() %>%
  arrange(desc(n)) %>%
  pull(sample_id) %>%
  head(5)

snvs_test <- snvs_tcga_brca %>%
  filter(sample_id %in% top_mutated_patients)
class(snvs_tcga_brca) <- c("cevo_snvs", class(snvs_test))


cna_hg19 <- read_tsv("/mnt/dane/data/brca_tcga_pan_can_atlas_2018/data_cna_hg19.seg")
cnvs_tcga_brca <- cna_hg19 |>
  transmute(
    sample_id = ID,
    chrom,
    start = loc.start,
    end = loc.end,
    log_ratio = NA_real_,
    BAF = NA_real_,
    total_cn = NA_real_,
    major_cn = NA_real_,
    minor_cn = NA_real_,
    seg_mean = seg.mean
  )

# usethis::use_data(cnvs_tcga_brca, overwrite = TRUE)
cnvs_test <- cnvs_tcga_brca %>%
  filter(sample_id %in% top_mutated_patients)

samples_data <- tibble(
  sample_id = unique(snvs_test$sample_id),
  patient_id = sample_id,
  sample = "tumor"
)

tcga_brca_test <- init_cevodata("TCGA BRCA test data", genome = "hg37") |>
  add_SNV_data(snvs_test, name = "TCGA") |>
  add_CNV_data(cnvs_test, data = "TCGA") |>
  add_sample_data(samples_data)
  # run_cevomod()

usethis::use_data(tcga_brca_test, overwrite = TRUE)
