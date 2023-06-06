
## --------------------------------- Data ------------------------------------

data("test_data", package = "cevoDatasets")
test_data

sample_data <- tibble(
  sample_id = test_data$metadata$sample_id,
  purity = 1
)

cnvs <- tribble(
  ~sample_id,  ~chrom, ~start, ~end,  ~total_cn, ~major_cn, ~minor_cn, ~normal_cn,
   "Sample 1",  "chr1", 1,      1000,  2,         1,         1,         2,
   "Sample 1",  "chr1", 1001,   2000,  4,         3,         1,         2,
   "Sample 1",  "chr1", 2001,   3000,  3,         2,         1,         2,
   "Sample 1",  "chr1", 3001,   4000,  2,         1,         1,         2,
   "Sample 2",  "chr1", 1,      1000,  2,         1,         1,         2,
   "Sample 2",  "chr1", 1001,   2000,  4,         3,         1,         2,
   "Sample 2",  "chr1", 2001,   3000,  3,         2,         1,         2,
   "Sample 2",  "chr1", 3001,   4000,  2,         1,         1,         2,
   "Sample 3",  "chr1", 1,      2000,  2,         1,         1,         2,
   "Sample 3",  "chr1", 2001,   4000,  4,         3,         1,         2,
   "Sample 4",  "chr1", 1,      4000,  2,         1,         1,         2
)



test_data <- test_data |>
  add_CNV_data(cnvs) |>
  add_sample_data(sample_data)

use_data(test_data, overwrite = TRUE)


## ------------------------ Fitted -------------------------------------------

withr::with_seed(
  123,
  test_data_fitted <- test_data |>
    prepare_SNVs() |>
    fit_powerlaw_tail_fixed() |>
    fit_subclones() |>
    fit_powerlaw_tail_optim() |>
    fit_subclones()
)

use_data(test_data_fitted, overwrite = TRUE)
