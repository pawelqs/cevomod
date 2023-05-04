

#' Fit subclonal distributions to neutral model residuals
#'
#' @param object object
#' @param N numbers of clones to for models
#' @param powerlaw_model_name residual of which powerlaw model to use?
#'   powerlaw_fixed/powerlaw_optim
#' @param snvs_name which snvs to to use?
#' @param method clustering method to use: BMix and mclust are currently supported.
#'   While mclust is a 3-4 times faster method, the BMix method is more accurate
#'   and usually fast enough.
#' @param upper_VAF_limit ignore variants with f higher than
#' @param verbose verbose?
#' @name fit_subclones
NULL



#' @rdname fit_subclones
#' @export
fit_subclones <- function(object,
                          N = 1:3,
                          powerlaw_model_name = active_models(object),
                          snvs_name = default_SNVs(object),
                          method = "BMix",
                          upper_VAF_limit = 0.75,
                          verbose = TRUE) {
  if (method == "BMix") {
    object <- object |>
      fit_subclones_bmix(N, powerlaw_model_name, snvs_name, upper_VAF_limit, verbose)
  } else if (method == "mclust") {
    object <- object |>
      fit_subclones_mclust(N, powerlaw_model_name, snvs_name, upper_VAF_limit, verbose)
  }

  object
}
