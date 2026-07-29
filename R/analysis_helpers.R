#' Collapse technical replicates
#'
#' @keywords internal

collapse_technical_replicates <- function(
    raw_data,
    technical_replicates = c("average", "error"),
    technical_sd_threshold = 0.2
) {
  technical_replicates <- match.arg(technical_replicates)

  prepared_data <- raw_data |>
    dplyr::transmute(
      Sample,
      Target,
      Cq = suppressWarnings(as.numeric(Cq)),
      `Amp Status`,
      valid_cq = dplyr::if_else(
        `Amp Status` == "Amp",
        Cq,
        NA_real_
      )
    )

  replicate_counts <- prepared_data |>
    dplyr::count(
      Sample,
      Target,
      name = "technical_n"
    )
  duplicated_pairs <- replicate_counts |>
    dplyr::filter(technical_n > 1L)

  if (nrow(duplicated_pairs) == 0L) {
    cli::cli_inform("No technical replicates were provided.")
  } else if (technical_replicates == "error") {
    cli::cli_abort(
      c(
        "Technical replicates were detected.",
        "x Repeated Sample/Target pairs: {nrow(duplicated_pairs)}"
      )
    )
  } else {
    cli::cli_inform(
      "Averaging technical replicates for {nrow(duplicated_pairs)} Sample/Target pair{?s}."
    )
  }

  collapsed_data <- prepared_data |>
    dplyr::group_by(
      Sample,
      Target
    ) |>
    dplyr::summarise(
      technical_n = dplyr::n(),
      technical_valid_n = sum(!is.na(valid_cq)),
      technical_not_determined_n = sum(is.na(valid_cq)),
      technical_sd = dplyr::if_else(
        technical_valid_n > 1L,
        stats::sd(valid_cq, na.rm = TRUE),
        NA_real_
      ),
      Cq = dplyr::if_else(
        technical_valid_n > 0L,
        mean(valid_cq, na.rm = TRUE),
        NA_real_
      ),
      `Amp Status` = dplyr::if_else(
        technical_valid_n > 0L,
        "Amp",
        "No Amp"
      ),
      .groups = "drop"
    )

  high_sd <- collapsed_data |>
    dplyr::filter(
      !is.na(technical_sd),
      technical_sd > technical_sd_threshold
    )
  if (nrow(high_sd) > 0L) {
    affected_pairs <- paste0(
      high_sd$Sample,
      " / ",
      high_sd$Target,
      " (SD = ",
      format(round(high_sd$technical_sd, 3), nsmall = 3),
      ")"
    )
    cli::cli_warn(
      c(
        "Technical-replicate SD exceeded {technical_sd_threshold} Cq.",
        "!" = "{affected_pairs}"
      )
    )
  }

  collapsed_data
}

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
