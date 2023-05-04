
data("test_data", package = "cevoDatasets")

test_data
use_data(test_data, overwrite = TRUE)


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
