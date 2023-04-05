
#' Fit Mobster models
#' @param object object
#' @param ... other args passed to mobster::fit_mobster()
#' @export
fit_mobster <- function(object, ...) {
  UseMethod("fit_mobster")
}


#' @rdname fit_mobster
#' @export
fit_mobster.cevodata <- function(object, ...) {
  rlang::check_installed("mobster")
  SNVs(object) |>
    select("sample_id", "VAF") |>
    nest_by(.data$sample_id) |>
    deframe() |>
    map(safely(mobster::mobster_fit), ...)
}
