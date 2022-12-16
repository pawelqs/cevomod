
# Get selected mutations based on 2 sample models
get_selected_mutations <- function(object, ...) {
  # As for now for one patient only
  if (!were_subclonal_models_fitted(object)) {
    stop("Fit subclonal models first!")
  }

  samples_data <- object$metadata |>
    select(sample_id:sample)
  rowsample <- samples_data$sample[[1]]
  colsample <- samples_data$sample[[2]]

  mutations_mat <- object |>
    get_SNVs_2d_matrix(rows_sample = rowsample, cols_sample = colsample, bins = 100)

  intervals <- tibble(
    interval = rownames(mutations_mat),
    centers = get_interval_centers(interval)
  )

  VAF_zero_pred <- tibble(
    VAF = 0, binom_pred = 0,
    sample = samples_data$sample
  )
  binom_predictions <- get_residuals(object) |>
    select(sample_id, VAF, binom_pred) |>
    left_join(samples_data, by = "sample_id") |>
    select(-sample_id)
  binom_predictions <- VAF_zero_pred |>
    bind_rows(binom_predictions) |>
    pivot_wider(names_from = "sample", values_from = "binom_pred") |>
    select(all_of(c("VAF", rowsample, colsample))) |>
    set_names(c("VAF", "rowsample", "colsample"))

  predictions_by_interval <- binom_predictions |>
    rebinarize_distribution(VAFs = intervals$centers) |>
    bind_cols(intervals)
  row_predictions <- deframe(predictions_by_interval[, c("interval", "rowsample")])
  col_predictions <- deframe(predictions_by_interval[, c("interval", "colsample")])

  limits <- init_MC_simulation_limits(mutations_mat, row_predictions, col_predictions)
  limit_ranges(limits) |> sum()

  # limits$upper[1:5, 1:5]
  mc_arr <- run_MC_simulation(upper_limits = limits$upper, iters = 2000)
  # mc_arr[, , 1]
  ev_res <- evaluate_MC_runs(mc_arr, row_predictions, col_predictions)
  plot(ev_res)

  top_models <- ev_res$err |>
    filter(rank < 0.01) |>
    arrange(err)
  # selected_solutions <- extract_models_to_tibble(mc_arr, top_models$i)
  plot_predictions_vs_fits(row_predictions, ev_res$rsums[top_models$i, ])
  plot_predictions_vs_fits(col_predictions, ev_res$csums[top_models$i, ])

  # average
  mean_top_solution <- average_solutions(mc_arr, top_models$i)
  plot_predictions_vs_fits(row_predictions, top_solution_rsums)
  plot_predictions_vs_fits(col_predictions, top_solution_csums)
  top_model_ev <- evaluate_MC_runs(mean_top_solution, row_predictions, col_predictions)
  top_model_ev
  plot(mean_top_solution)

  # Another round
  new_limits <- tune_limits(limits, mc_arr[, , top_models$i], row_predictions, col_predictions)
  limit_ranges(limits) |> sum()
  limit_ranges(new_limits) |> sum()

  # Loop
  limits <- init_MC_simulation_limits(mutations_mat, row_predictions, col_predictions)
  limit_ranges(limits) |> sum()
  for (i in 1:5) {
    limit_ranges(limits) |> sum() |> print()
    mc_arr <- run_MC_simulation(upper_limits = limits$upper, lower_limits = limits$lower, iters = 20000)
    ev_res <- evaluate_MC_runs(mc_arr, row_predictions, col_predictions)

    top_models <- ev_res$err |>
      filter(rank < 0.01) |>
      arrange(err)

    mean_top_solution <- average_solutions(mc_arr, top_models$i)
    top_model_ev <- evaluate_MC_runs(mean_top_solution, row_predictions, col_predictions)
    top_model_ev |> print()

    limits <- tune_limits(limits, mc_arr[, , top_models$i], row_predictions, col_predictions)
  }
}


init_MC_simulation_limits <- function(mutations_mat, row_predictions, col_predictions) {
  upper_limits <- mutations_mat
  zero_interval <- rownames(mutations_mat)[[1]]
  for (row in rownames(mutations_mat)) {
    for (col in colnames(mutations_mat)) {
      if (row == zero_interval && col == zero_interval) {
        upper_limits[row, col] <- 0
      } else if (row == zero_interval) {
        upper_limits[row, col] <- min(mutations_mat[row, col], col_predictions[col] * 1.1)
      } else if (col == zero_interval) {
        upper_limits[row, col] <- min(mutations_mat[row, col], row_predictions[row] * 1.1)
      } else {
        upper_limits[row, col] <- min(mutations_mat[row, col], row_predictions[row] * 1.1, col_predictions[col] * 1.1)
      }
    }
  }
  upper_limits <- round(upper_limits)
  lower_limits <- upper_limits
  lower_limits[, ] <- 0
  limits <- list(lower = lower_limits, upper = upper_limits)
  class(limits) <- c("cevo_MC_sim_limits", class(limits))
  limits
}


limit_ranges <- function(limits, zero.rm = TRUE) {
  zero_limits <- limits$upper == 0
  limits_diff <- limits$upper - limits$lower
  c(limits_diff[!zero_limits])
}


tune_limits <- function(limits, mc_arr, row_predictions, col_predictions) {
  upper_limits <- limits$upper
  lower_limits <- limits$lower
  for (row in rownames(upper_limits)) {
    for (col in colnames(upper_limits)) {
      max_val <- max(mc_arr[row, col, ])
      min_val <- min(mc_arr[row, col, ])
      # mean_val <- mean(mc_arr[row, col, ])
      # sd <- sd(mc_arr[row, col, ])
      # max_val <- mean_val + 2 * sd
      # min_val <- mean_val - 2 * sd
      # upper
      if (max_val < row_predictions[row] || max_val < col_predictions[col]) {
        upper_limits[row, col] <- upper_limits[row, col]
      } else {
        upper_limits[row, col] <- max_val
      }
      # lower
      if (min_val > row_predictions[row] || min_val > col_predictions[col]) {
        # lower_limits[row, col] <- 0
        lower_limits[row, col] <- lower_limits[row, col]
      } else {
        lower_limits[row, col] <- min_val
      }
    }
  }
  upper_limits <- round(upper_limits)
  lower_limits <- round(lower_limits)
  limits <- list(lower = lower_limits, upper = upper_limits)
  class(limits) <- c("cevo_MC_sim_limits", class(limits))
  limits
}


run_MC_simulation <- function(upper_limits, lower_limits = NULL, iters = 2000) {
  mc_arr <- array(
    dim = c(nrow(upper_limits), ncol(upper_limits), iters),
    dimnames = list(rownames(upper_limits), colnames(upper_limits), 1:iters)
  )
  for (row in rownames(upper_limits)) {
    for (col in colnames(upper_limits)) {
      mc_arr[row, col, ] <- runif(
        n = iters,
        min = if (is.null(lower_limits)) 0 else lower_limits[row, col],
        max = upper_limits[row, col]
      )
    }
  }
  mc_arr <- round(mc_arr)
  mc_arr
}


extract_models_to_tibble <- function(mc_arr, which) {
  which |>
    set_names(which) |>
    map(~mc_arr[, , .x]) |>
    map(as.data.frame) |>
    map(rownames_to_column, "rowsample") |>
    map(~pivot_longer(.x, -rowsample, names_to = "colsample", values_to = "N")) |>
    bind_rows(.id = "i")
}


average_solutions <- function(mc_arr, which) {
  mean_solution <- array(
    dim = dim(mc_arr)[1:2],
    dimnames = dimnames(mc_arr)[1:2]
  )
  for (row in rownames(mean_solution)) {
    for (col in colnames(mean_solution)) {
      mean_solution[row, col] <- mean(mc_arr[row, col, which])
    }
  }
  mean_solution
  class(mean_solution) <- c("non_neutral_2d_fit", class(mean_solution))
  mean_solution
}


# ----------------------- Models evaluation -----------------------------------

#' @export
evaluate_MC_runs <- function(mc_arr, rowsums_pred, colsums_pred) {
  UseMethod("evaluate_MC_runs")
}


#' @export
evaluate_MC_runs.array <- function(mc_arr, rowsums_pred, colsums_pred) {
  iter <- dim(mc_arr)[[3]]
  err <- rep(NA_real_, iter)
  rsums <- matrix(NA, nrow = iter, ncol = dim(mc_arr)[[1]])
  csums <- matrix(NA, nrow = iter, ncol = dim(mc_arr)[[2]])
  # non_zero_sums <- matrix(NA, nrow = iter, ncol = length(row_predictions) - 1)
  for (i in 1:iter) {
    rsums[i, ] <- rowSums(mc_arr[, , i])
    csums[i, ] <- colSums(mc_arr[, , i])
    non_zero_sums <- c(rsums[i, -1], csums[i, -1])
    non_zero_preds <- c(rowsums_pred[-1], colsums_pred[-1])
    err[i] <- sum((non_zero_sums - non_zero_preds)^2)
  }
  err <- tibble(err, i = 1:iter) |>
    mutate(rank = percent_rank(err))
  res <- lst(err, rsums, csums)
  class(res) <- c("cevo_MC_solutions_eval", class(res))
  res
}


#' @export
evaluate_MC_runs.matrix <- function(mc_arr, rowsums_pred, colsums_pred) {
  rsums <- rowSums(mc_arr)
  csums <- colSums(mc_arr)
  non_zero_sums <- c(rsums[-1], csums[-1])
  non_zero_preds <- c(rowsums_pred[-1], colsums_pred[-1])
  err <- sum((non_zero_sums - non_zero_preds)^2)
  rsq_rows <- rsq_vec(rowsums_pred[-1], rsums[-1])
  rsq_cols <- rsq_vec(colsums_pred[-1], csums[-1])
  rsq_total <- rsq_vec(non_zero_preds, non_zero_sums)
  res <- lst(err, rsq_rows, rsq_cols, rsq_total, rsums, csums)
  class(res) <- c("non_neutral_2d_fit_eval", class(res))
  res
}


#' @export
print.non_neutral_2d_fit_eval <- function(x, ...) {
  cli::cat_line("MSE:     ", round(x$err, digits = 2))
  cli::cat_line("Rows R^2:    ", round(x$rsq_rows, digits = 2))
  cli::cat_line("Cols R^2:    ", round(x$rsq_cols, digits = 2))
  cli::cat_line("Total R^2:   ", round(x$rsq_total, digits = 2))
  cli::cat_line("Row sums:  ", paste0(x$rsums[1:5], collapse = ", "), "...")
  cli::cat_line("Col sums:  ", paste0(x$csums[1:5], collapse = ", "), "...")
}


#' @export
print.cevo_MC_solutions_eval <- function(x, ...) {
  MSE_range_str <- x$err$err |>
    range() |>
    round(digits = 0) |>
    paste0(collapse = " - ")
  cli::cat_line("Iters:     ", nrow(x$err))
  cli::cat_line("MSE range: ", MSE_range_str)
}


#' @export
plot.cevo_MC_solutions_eval <- function(x, ...) {
  hist(x$err$err)
}


# -------------------------------- Plots --------------------------------------

plot_predictions_vs_fits <- function(predictions, fits) {
  plot(predictions)
  if (is.null(dim(fits))) {
    points(fits[-1], col = "red")
  } else {
    for (i in 1:nrow(fits)) {
      points(fits[i, -1], col = "red")
    }
  }
}


#' @export
plot.non_neutral_2d_fit <- function(x, ...) {
  pheatmap::pheatmap(
    x, cluster_rows = FALSE, cluster_cols = FALSE
  )
}

# join_patient_models <- function(object, ...) {
#   if (!were_subclonal_models_fitted(object)) {
#     stop("Fit subclonal models first!")
#   }
#
#   binom_models <- get_models(object) |>
#     left_join(object$metadata, by = "sample_id")
#   total_list_of_clones <- binom_models |>
#     filter(component != "Neutral tail") |>
#     select(patient_id, sample, component, N_mutations, cellularity) |>
#     mutate(clone_id = str_c(component, sample, sep = "_"), .before = "sample")
#
#   snvs <- SNVs(object) |>
#     select(sample_id, chrom:alt, VAF, DP, impact) |>
#     left_join(object$metadata |> select(patient_id:sample), by = "sample_id") |>
#     unite(mutation_id, chrom:alt, sep = "-")
#
#   mut_predictions <- snvs |>
#     classify_SNVs(get_residuals(object)) |>
#     select(patient_id, mutation_id, sample, starts_with("Clone"), starts_with("Subclone"))
#
#   mut_predictions_wide <- mut_predictions |>
#     pivot_wider(id_cols = patient_id:mutation_id, names_from = sample, values_from = Clone:`Subclone 2`) |>
#     binarise_mut_predictions()
#
#   clone_links <- mut_predictions |>
#     select(-mutation_id) |>
#     nest_by(patient_id) |>
#     deframe() |>
#     map(count_clone_linking_mutations) |>
#     bind_rows(.id = "patient_id") |>
#     separate(cloneid1, into = c("clone1", "sample1"), sep = "_", remove = FALSE) |>
#     separate(cloneid2, into = c("clone2", "sample2"), sep = "_", remove = FALSE) |>
#     filter(sample1 != sample2)
#
#   x <- clone_links |>
#     left_join(total_list_of_clones |> select(patient_id, cloneid1 = clone_id, N_mutations1 = N_mutations)) |>
#     left_join(total_list_of_clones |> select(patient_id, cloneid2 = clone_id, N_mutations2 = N_mutations)) |>
#     mutate(
#       clone1_link_frct = common_mutations / N_mutations1,
#       clone2_link_frct = common_mutations / N_mutations2,
#       link_type = case_when(
#         clone1_link_frct > 0.8 & clone2_link_frct > 0.8 ~ "strong link",
#         clone1_link_frct > 0.8 & clone2_link_frct < 0.8 ~ "clone1 part of clone2",
#         clone1_link_frct < 0.8 & clone2_link_frct > 0.8 ~ "clone2 part of clone1",
#         TRUE ~ "no link"
#       )
#     )
#
#   x <- total_list_of_clones |>
#     left_join(clone_links, by = c("patient_id", "sample", "clone"))
#
#   clones <- tibble(
#     patient_id,
#     clone_id,
#     sample_id,
#     clone_name,
#     cellularity,
#     N_mutations
#   )
#
#   object$models[["clones"]] <- clones
#   # re-classify mutations
#   # draw tree
#   # get evo params
# }
#
#
#
# binarise_mut_predictions <- function(mut_predictions) {
#   probability_cols <- colnames(mut_predictions) |>
#     str_subset(pattern = "[Cc]lone")
#   mut_predictions |>
#     group_by(.data$patient_id) |>
#     mutate_at(probability_cols, replace_na, 0) |>
#     mutate_at(probability_cols, ~.x/max(.x, na.rm = TRUE)) |>
#     mutate_at(probability_cols, ~.x >= 0.0001) |>
#     ungroup()
# }
#
#
# # # binary_links_df <- mut_predictions |>
# # #   filter(patient_id == "AMLRO-10") |>
# # #   select(-patient_id, -mutation_id)
# #
# # count_clone_linking_mutations <- function(binary_links_df) {
# #   variables <- colnames(binary_links_df)
# #   res <- utils::combn(variables, 2, simplify = FALSE) |>
# #     map(set_names, c("cloneid1", "cloneid2")) |>
# #     bind_rows() |>
# #     mutate(
# #       common_mutations = map2_int(
# #         .data$cloneid1, .data$cloneid2,
# #         ~sum(binary_links_df[[.x]] & binary_links_df[[.y]])
# #       )
# #     ) |>
# #     filter(!is.na(common_mutations))
# #   res
# # }
# # #   object$models$binomial_models |>
# # #     left_join(object$metadata) |>
# # #     transmute(
# # #       patient_id,
# # #       full_clone_id_1 = str_c(sample, clone, sep = "_")
# # #     ) |>
# # #     nest_by(patient_id) |>
# # #     deframe() |>
# # #     map(~expand_grid(.x, full_clone_id_2 = .x$full_clone_id_1))
# # #   expand_grid(full_clone_id_2 = full_clone_id)
# # #
# # #   connections <- object$metadata |>
# # #     select(patient_id, sample) |>
# # #     nest_by(patient_id) |>
# # #     deframe() |>
# # #     map("sample") |>
# # #     map(sort) |>
# # #     map(function(samples) {
# # #       samples |>
# # #         combn(2) |>
# # #         t() |>
# # #         as.data.frame() |>
# # #         as_tibble() |>
# # #         set_names(c("sample1", "sample2"))
# # #     }) |>
# # #     bind_rows(.id = "patient_id")
# # #   expand_grid(sample2 = sample)
# # #
# # #
# # #   corr <- classifications |>
# # #     select(-Neutral) |>
# # #     pivot_wider(id_cols = patient_id:mutation_id, names_from = sample, values_from = Clone:`Subclone 2`) |>
# # #     select(-mutation_id) |>
# # #     nest() |>
# # #     deframe() |>
# # #     map(corrr::correlate, quiet = TRUE) |>
# # #     map(corrr::shave) |>
# # #     map(corrr::stretch) |>
# # #     map(~separate(.x, into = ""))
# # #   map(~filter(.x, !is.na(.data$r)))
# # #   map(column_to_rownames, "term") |>
# # #     map(janitor::remove_empty, which = c("rows", "cols")) |>
# # #     map(rownames_to_column, "term")
# # #   select(-data)
# # #   corr
# # #
# # #   p <- corr$corr_res |>
# # #     map(function(x) {
# # #       x <- x |>
# # #         keep(~!all(is.na(.x)))
# # #       x[is.na(x)] <- 0
# # #       x |>
# # #         column_to_rownames("term")
# # #     }) |>
# # #     map(pheatmap::pheatmap)
# # # }
#
#
# classify_SNVs <- function(snvs, residuals) {
#   probabilities <- get_VAF_clonal_probabilities_tbl(residuals)
#   snvs |>
#     mutate(VAF_chr = as.character(round(.data$VAF, digits = 2))) |>
#     left_join(probabilities, by = c("sample_id", "VAF_chr")) |>
#     select(-.data$VAF_chr)
# }
#
#
# get_VAF_clonal_probabilities_tbl <- function(residuals) {
#   probabilities <- residuals |>
#     mutate(VAF = as.character(.data$VAF)) |>
#     select(
#       .data$sample_id, .data$VAF,
#       Neutral = .data$neutral_pred,
#       starts_with("Clone"),
#       starts_with("Subclone"),
#       .data$model_pred
#     ) |>
#     transmute(
#       .data$sample_id,
#       VAF_chr = .data$VAF,
#       across(c("Neutral", starts_with("Clone"), starts_with("Subclone")), ~.x/model_pred)
#     )
#   probabilities
# }
