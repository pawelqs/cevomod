# linear_models <- models
#
# fit_binomial_models <- function(sfs, linear_models) {
#   linear_models <- linear_models |>
#     slice(1)
#   err <- estimate_sampling_rate(sfs, linear_models) |>
#     left_join(linear_models) |>
#     filter(VAF >= to)
# }
#
#
# ggplot(err, aes(VAF, -err)) +
#   geom_point(size = 0.2) +
#   facet_wrap(~sample_id)
