## code to prepare `sysdata` goes here
library(tidyverse)


driver_genes <- openxlsx::read.xlsx(
    "data-raw/driver_genes.xlsx",
    sheet = "Table S1", startRow = 4
  ) |>
  as_tibble() |>
  rename(
    Type = `Tumor.suppressor.or.oncogene.prediction.(by.20/20+)`
  )


variant_classes <- openxlsx::getSheetNames("data-raw/variant_classification.xlsx")
variant_classification <- variant_classes %>%
  set_names(variant_classes) %>%
  map(~openxlsx::read.xlsx("data-raw/variant_classification.xlsx", sheet = .x, colNames = FALSE)) %>%
  map(set_names, "Variant_Classification") %>%
  bind_rows(.id = "variant_classification")


usethis::use_data(
  driver_genes, variant_classification,
  overwrite = TRUE, internal = TRUE
)
usethis::use_data(driver_genes, overwrite = TRUE)
usethis::use_data(variant_classification, overwrite = TRUE)
