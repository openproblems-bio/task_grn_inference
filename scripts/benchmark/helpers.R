
`%|%` <- function(x, y) {
  ifelse(is.na(x), y, x)
}
aggregate_scores <- function(score) {
  mean(pmin(1, pmax(0, score)) %|% 0)
}
label_time <- function(time) {
  case_when(
    time < 1e-5 ~ "0s",
    time < 1 ~ "<1s",
    time < 60 ~ paste0(floor(time), "s"),
    time < 3600 ~ paste0(floor(time / 60), "m"),
    time < 3600 * 24 ~ paste0(floor(time / 3600), "h"),
    time < 3600 * 24 * 7 ~ paste0(floor(time / 3600 / 24), "d"),
    !is.na(time) ~ ">7d",
    TRUE ~ NA_character_
  )
}
label_memory <- function(x, include_mb = FALSE) {
  case_when(
    x < 1e9 ~ "<1G",
    x < 1e12 ~ paste0(round(x / 1e9), "G"),
    !is.na(x) ~ ">1T",
    TRUE ~ NA_character_
  )
}
normalize_scores <- function(
  df,
  groups = c("dataset_id", "metric_id"),
  metric_value = "value",
  maximize = "maximize",
  rescale_only_with_controls = FALSE
) {
  df2 <- df
  if (rescale_only_with_controls) {
    df2 <- df2 %>% filter(is_baseline)
  }
  norm <- df2 %>%
    group_by(!!!syms(groups)) %>%
    summarise(
      min = min(!!sym(metric_value), na.rm = TRUE),
      max = max(!!sym(metric_value), na.rm = TRUE),
      .groups = "drop"
    )
  scaled <- df %>%
    left_join(norm, by = groups) %>%
    mutate(
      score = (!!sym(metric_value) - min) / (max - min),
      score = ifelse(!!sym(maximize), score, 1 - score),
      score = score %|% 0
    )

  list(norm = norm, scaled = scaled)
}