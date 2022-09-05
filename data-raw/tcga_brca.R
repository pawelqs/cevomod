## code to prepare `tcga_brca` dataset goes here

# Go here -> https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
# Download the data using the Download Icon next to the dataset name
# Unpack
# Fix the path below

library(tidyverse)

mutations <- read_tsv("/mnt/dane/data/brca_tcga_pan_can_atlas_2018/data_mutations.txt")
glimpse(mutations)


tcga_brca <- mutations %>%
  left_join(variant_classification) %>%
  transmute(
    sample_id = Tumor_Sample_Barcode,
    chrom = Chromosome,
    start = Start_Position,
    gene_symbol = Hugo_Symbol,
    ref = Reference_Allele,
    alt = Tumor_Seq_Allele2,
    ref_counts = t_ref_count,
    alt_counts = t_alt_count,
    VAF = t_alt_count / (t_alt_count + t_ref_count),
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
class(tcga_brca) <- c("cevo_SNVs_tbl", class(tcga_brca))

usethis::use_data(tcga_brca, overwrite = TRUE)

top_mutated_patients <- tcga_brca %>%
  group_by(sample_id) %>%
  count() %>%
  arrange(desc(n)) %>%
  pull(sample_id) %>%
  head(5)

tcga_brca_test <- tcga_brca %>%
  filter(sample_id %in% top_mutated_patients)
class(tcga_brca) <- c("cevo_SNVs_tbl", class(tcga_brca_test))

usethis::use_data(tcga_brca_test, overwrite = TRUE)

