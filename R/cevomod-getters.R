
#' Get model names
#' @param object object
#' @export
get_model_names <- function(object) {
  UseMethod("get_model_names")
}


#' @export
get_model_names.cevodata <- function(object) {
  names(object$models)
}


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
    stop("Slot ", which, " is empty! Fit apropriate model first")
  }
  if (best_only && !is.null(models[["best"]])) {
    filter(models, .data$best)
  } else {
    models
  }
}


get_powerlaw_models <- function(object,
                                which = active_models(object),
                                best_only = TRUE,
                                ...) {
  models <- get_models(object, which)
  if ("cevo_powerlaw_models" %not in% class(models)) {
    stop(
      which, " is not a powerlaw model, required to fit subclones.",
      "Use another model"
    )
  }
  models
}


active_models <- function(object, ...) {
  if (is.null(object$active_models) | length(object$models) == 0) {
    stop("No models has been fitted yet!")
  }
  object$active_models
}
