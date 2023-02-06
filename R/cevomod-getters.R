
#' Get models from the object
#' @param object object to get the models from
#' @param ... other arguments
#' @export
get_models <- function(object, ...) {
  UseMethod("get_models")
}


#' @describeIn get_models Get models from cevodata object
#' @param which `chr` which models to get
#' @param best_only `lgl` return only the best models?
#' @export
get_models.cevodata <- function(object,
                                which = active_models(object),
                                best_only = TRUE,
                                ...) {
  models <- object$models[[which]]
  if (is.null(models)) {
    stop("Slot ", name, " is empty! Fit apropriate model first")
  }
  if (best_only && !is.null(models$best)) {
    filter(models, .data$best)
  } else {
    models
  }
}


active_models <- function(object, ...) {
  if (is.null(object$active_models) | length(object$models) == 0) {
    stop("No models has been fitted yet!")
  }
  object$active_models
}
