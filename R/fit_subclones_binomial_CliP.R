
#' @describeIn fit_subclones Fit subclonal distributions to neutral model
#'   residuals using CliP - Clonal structure identification through penalizing
#'   pairwise differences [(github)](https://github.com/wwylab/CliP). Requires
#'   that the Apptainer installed and
#'   - the path to the CliP.sif container image is provided,
#'   - or CliP.sif exists in the containers_dir ([set_containers_dir()])
#'   - or CliP.sif exists in the current working directory.
#'   The container should be build with [build_clip_container()].
#' @param clip_sif Apptainer image file with CliP. If NULL, cevomod seeks for
#'   the CliP.sif file in the containers_dir (if set before) and in the current
#'   working directory. The container should be build with [build_clip_container()].
#' @param clip_intput Path to store the CliP input files
#' @param clip_output Path to store the CliP output files
#' @export
fit_subclones_clip <- function(object,
                               powerlaw_model_name = active_models(object),
                               snvs_name = default_SNVs(object),
                               cnvs_name = default_CNVs(object),
                               upper_f_limit = 0.75,
                               clip_sif = NULL,
                               clip_input = file.path(tempdir(), "clip_input"),
                               clip_output = file.path(tempdir(), "clip_output"),
                               verbose = get_cevomod_verbosity()) {
  msg("Fitting binomial models using CllP", verbose = verbose)
  rlang::check_installed("readthis", reason = "to read CliP results")
  if (!is_apptainer_installed()) {
    stop("Apptainer needs to be installed to run fit subclones using CliP")
  }

  if (is.null(clip_sif)) {
    clip_sif <- find_clip_container()
  }

  powerlaw_models <- get_powerlaw_models(object, powerlaw_model_name)
  residuals <- get_residuals(object, models_name = powerlaw_model_name) |>
    filter(.data$f >= 0)
  sequencing_depth <- SNVs(object, which = snvs_name) |>
    get_local_sequencing_depths() |>
    transmute(.data$sample_id, .data$f, sequencing_DP = .data$median_DP)

  non_neutral_tail_mut_counts <- residuals |>
    mutate(n = if_else(.data$f > upper_f_limit, 0, round(.data$powerlaw_resid_clones))) |>
    select("sample_id", "f_interval", "n")
  snvs_to_cluster <- SNVs(object, which = snvs_name) |>
    nest_by(.data$sample_id, .data$f_interval) |>
    inner_join(non_neutral_tail_mut_counts, by = c("sample_id", "f_interval")) |>
    reframe(
      slice_sample(.data$data, n = .data$n)
    )

  clip_files <- object |>
    add_SNV_data(snvs_to_cluster, "snvs_to_cluster") |>
    to_clip(out_dir = clip_input, snvs_name = "snvs_to_cluster", cnvs_name = cnvs_name)

  pmap(
    clip_files,
    run_clip,
    clip_sif = clip_sif,
    out_dir = clip_output,
    verbose = verbose,
    .progress = TRUE
  )

  clip_res <- readthis::read_clip_best_lambda(clip_output)

  models <- clip_res$subclonal_structure |>
    transmute(
      sample_id,
      model = "subclones CliP",
      component = if_else(cluster_index == 0, "Clone", str_c("Subclone ", cluster_index)),
      N_mutations = num_SNV,
      cellularity = .data$cellular_prevalence/2,
      f = round(cellularity, digits = 2),
      lambda,
      best = best_lambda
    ) |>
    left_join(sequencing_depth, by = c("sample_id", "f")) |>
    select(-"f")

  best_models <- models |>
    filter(.data$best) |>
    nest_by(.data$sample_id, .key = "clones")

  clonal_predictions <- residuals |>
    select("sample_id", "f_interval", "f") |>
    nest_by(.data$sample_id, .key = "intervals") |>
    inner_join(best_models, by = "sample_id") |>
    reframe(get_binomial_predictions(.data$clones, .data$intervals)) |>
    select(-"f")

  residuals <- residuals |>
    select(-"model_resid") |>
    left_join(clonal_predictions, by = c("sample_id", "f_interval")) |>
    mutate(
      model_pred = .data$powerlaw_pred + .data$binom_pred,
      model_resid = .data$SFS - .data$model_pred
    )

  models <- powerlaw_models |>
    bind_rows(models) |>
    arrange(.data$sample_id, .data$best, .data$model)

  models_name <- paste0(powerlaw_model_name, "_subclones")
  resid_name <- paste0("residuals_", models_name)
  object$models[[models_name]] <- models
  object$misc[[resid_name]] <- residuals
  object$active_models <- models_name
  object
}



find_clip_container <- function() {
  containers_dir <- get_containers_dir()
  clip_sif <- if (!is.null(containers_dir)) {
    file.path(containers_dir, "CliP.sif") |>
      str_replace("//", "/")
  }
  if (file.exists(clip_sif)) {
    return(clip_sif)
  }

  clip_sif <- "CliP.sif"
  if (file.exists(clip_sif)) {
    return(clip_sif)
  }

  stop("Provide clip_sif path or use build_clip_container() to build it.")
}



run_clip <- function(sample_id, snvs_file, cnvs_file, purity_file, out_dir,
                     clip_sif = NULL,
                     clip_args = c("-O", out_dir, "--sample_id", sample_id),
                     verbose = TRUE) {
  rlang::check_installed("processx", reason = "to run CliP")
  if (!is_apptainer_installed()) {
    stop("Apptainer needs to be installed to run fit subclones using CliP")
  }

  args <- c(
    "exec", clip_sif,
    "python", "/opt/CliP/run_clip_main.py",
    snvs_file, cnvs_file, purity_file, clip_args
  )
  out <- processx::run("apptainer", args, echo = verbose, echo_cmd = verbose)

  invisible(out)
}

