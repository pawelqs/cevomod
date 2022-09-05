'%not in%' <- function(x,y)!('%in%'(x,y))


join_aes <- function(aes_default, aes_2) {
  aes <- c(as.list(aes_default), aes_2[!names(aes_2) %in% names(aes_default)])
  class(aes) <- 'uneval'
  aes
}


complete_missing_VAF_levels <- function(dt, fill, digits = 2) {
  VAF <- NULL
  group_variables <- group_vars(dt)
  dt %>%
    mutate(VAF = round(VAF, digits = digits)) %>%
    mutate(
      VAF = VAF %>%
        as.character() %>%
        parse_factor(levels = as.character(seq(0, 1, by = 1/(10^digits))))
    ) %>%
    complete(VAF, fill = fill) %>%
    mutate(
      VAF = VAF %>%
        as.character() %>%
        parse_double()
    ) %>%
    group_by(!!!syms(group_variables))
}
