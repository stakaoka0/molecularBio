#' Analyze qPCR data
#'
#' Reads QuantStudio results and sample metadata, then calculates relative
#' expression using a delta-Cq workflow and geometric-mean normalization
#' across the specified reference genes.
#'
#' @param path_to_qPCR_res Path to a QuantStudio CSV file. Required columns are
#'   `Sample`, `Target`, `Cq`, and `Amp Status`.
#' @param path_to_group_label Path to a sample metadata CSV file. Required
#'   columns are `Sample`, `Group`, and `Group_order`. An optional
#'   `color_scheme` column controls grouped plot colors.
#' @param controls Character vector of reference-gene names.
#' @param technical_replicates How repeated `Sample`-`Target` combinations are
#'   handled: `"average"` averages valid Cq values, while `"error"` rejects
#'   repeated combinations.
#' @param technical_sd_threshold Non-negative Cq SD threshold used to warn
#'   about variable technical replicates.
#' @param out_dir_path Output directory used when `save_results = TRUE`.
#' @param save_results Logical; save the analyzed table as
#'   `analyzed_qPCR_res.tsv`.
#'
#' @return A `qpcr_analysis` object.
#' @export

analyze_qPCR <- function(
  path_to_qPCR_res,
  path_to_group_label,
  controls = c("18S", "GAPDH", "RPL13A"),
  technical_replicates = c("average", "error"),
  technical_sd_threshold = 0.2,
  out_dir_path = NULL,
  save_results = FALSE
) {
  technical_replicates <- match.arg(technical_replicates)
  if (!is.numeric(technical_sd_threshold) ||
      length(technical_sd_threshold) != 1L ||
      is.na(technical_sd_threshold) ||
      !is.finite(technical_sd_threshold) ||
      technical_sd_threshold < 0) {
    cli::cli_abort(
      "{.arg technical_sd_threshold} must be one non-negative finite number."
    )
  }
  if (!is.logical(save_results) || length(save_results) != 1L ||
      is.na(save_results)) {
    cli::cli_abort("{.arg save_results} must be `TRUE` or `FALSE`.")
  }
  if (save_results && is.null(out_dir_path)) {
    cli::cli_abort(
      "{.arg out_dir_path} is required when {.arg save_results} is `TRUE`."
    )
  }

  raw_data <- readr::read_csv(
    path_to_qPCR_res,
    col_types = readr::cols()
  )

  group_label <- readr::read_csv(
    path_to_group_label,
    col_types = readr::cols()
  )
  if (!"color_scheme" %in% names(group_label)) {
    group_label$color_scheme <- 1L
  }

  validate_qpcr_inputs(
    raw_data,
    group_label,
    controls
  )
  raw_data <- collapse_technical_replicates(
    raw_data,
    technical_replicates = technical_replicates,
    technical_sd_threshold = technical_sd_threshold
  )

  base_group <- group_label |>
    dplyr::filter(Group_order == 1) |>
    dplyr::pull(Group) |>
    unique()

  summary_data <- raw_data |>
    dplyr::left_join(
      group_label |>
        dplyr::select(
          Sample,
          Group,
          Group_order
      ),
      by = "Sample"
    )

  # ΔCt
  # RQ
  # normalization factor
  # normalized expression
  final_data <- build_expression_table(
    summary_data,
    controls,
    base_group
  )

  group_levels <- group_label %>%
    dplyr::distinct(
      Group,
      Group_order
    ) %>%
    dplyr::arrange(
      Group_order
    ) %>%
    dplyr::pull(Group)

  if (save_results) {

    dir.create(
      out_dir_path,
      recursive = TRUE,
      showWarnings = FALSE
    )

    readr::write_tsv(
      final_data,
      file.path(
        out_dir_path,
        "analyzed_qPCR_res.tsv"
      )
    )
  }

  return(
    structure(
      list(
        data = final_data,
        controls = controls,
        reference_group = base_group,
        group_levels = group_levels,

        group_info = group_label %>%
          dplyr::distinct(
            Group,
            Group_order,
            color_scheme
          )
      ),
      class = "qpcr_analysis"

    )
  )
}
