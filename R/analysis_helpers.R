#' Build expression table
#'
#' @keywords internal

build_expression_table <- function(
    summary_data,
    controls,
    base_group
) {
  # Mean Ct per group
  summarize_data <- summary_data |>
    dplyr::group_by(
      Group,
      Target
    ) |>
    dplyr::mutate(
      Cq_mean = mean(
        Cq,
        na.rm = TRUE
      )
    )
  # Reference group means
  base_mean <- summarize_data |>
    dplyr::filter(
      Group == base_group
    ) |>
    dplyr::distinct(
      Target,
      Cq_mean
    )

  # ΔCt and RQ

  summarize_data_dct <- summarize_data |>
    dplyr::left_join(
      base_mean,
      by = "Target",
      suffix = c("", "_base")
    ) |>
    dplyr::mutate(
      dCt = Cq_mean_base - Cq,
      RQ = 2^dCt
    )

  # Normalization factors

  norm_factors <- summarize_data_dct |>
    dplyr::filter(
      Target %in% controls
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      Sample
    ) |>
    dplyr::summarise(
      norm_factor = geo_mean(RQ),
      .groups = "drop"
    )

  # Final table

  summarize_data_dct |>
    dplyr::left_join(
      norm_factors,
      by = "Sample"
    ) |>
    dplyr::mutate(
      norm_exp = RQ / norm_factor,
      log2_norm_exp = log2(norm_exp)
    ) |>
    dplyr::group_by(
      Group,
      Target
    ) |>
    dplyr::mutate(
      group_exp = geo_mean(norm_exp),
      log2_group_exp = log2(group_exp),
      not_determined_n = sum(
        is.na(Cq)
      ),
      n = sum(
        !is.na(Cq)
      ),
      total_n = not_determined_n + n,
      sd_log2 = dplyr::if_else(
        n > 1,
        stats::sd(
          log2_norm_exp,
          na.rm = TRUE
        ),
        NA_real_
      ),
      sem = sd_log2 / sqrt(n),
      upper_error =
        2^(log2_group_exp + sem),
      lower_error =
        2^(log2_group_exp - sem)
    ) |>
    dplyr::ungroup()

  group_order_lookup <- summary_data |>
    dplyr::distinct(
      Group,
      Group_order
    )

  final_data <- summarize_data_dct %>%
    dplyr::left_join(
      norm_factors,
      by = "Sample"
    ) %>%
    dplyr::mutate(
      norm_exp = RQ / norm_factor,
      log2_norm_exp = log2(norm_exp)
    ) %>%
    dplyr::group_by(
      Group,
      Group_order,
      Target
    ) %>%
    dplyr::mutate(
      group_exp = geo_mean(norm_exp),
      log2_group_exp = log2(group_exp),
      not_determined_n = sum(is.na(Cq)),
      n = sum(!is.na(Cq)),
      total_n = not_determined_n + n,
      sd_log2 = dplyr::if_else(
        n > 1,
        stats::sd(
          log2_norm_exp,
          na.rm = TRUE
        ),
        NA_real_
      ),
      sem = sd_log2 / sqrt(n),
      upper_error =
        2^(log2_group_exp + sem),
      lower_error =
        2^(log2_group_exp - sem)
    ) %>%
    dplyr::ungroup()

  return(final_data)
}
