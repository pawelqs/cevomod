suppressPackageStartupMessages(library(cevomod))
suppressPackageStartupMessages(library(cevoDatasets))
suppressPackageStartupMessages(library(furrr))
plan(multisession, workers = 4)

cli::cat_line("Running cevomod on cevoDatasets", col = "green4")

data <- list(
  AMLRO = AMLRO,
  OPUS_BRCA = OPUS_BRCA,
  OPUS_Larynx = OPUS_Larynx
)
data <- future_map(data, run_cevomod)

dir.create("~/.cevomod", showWarnings = FALSE)
readr::write_rds(data, "~/.cevomod/data.Rds")

cli::cat_line("cevomod results saved to ~/.cevomod/data.Rds", col = "green4")
