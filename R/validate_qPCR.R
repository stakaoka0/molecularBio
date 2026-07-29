# validate_qPCR.R
#' Validate qPCR inputs
#'
#' @keywords internal
#'
validate_qpcr_inputs <- function(
  raw_data,
  group_label,
  controls
) {
  required_qpcr_cols <- c(
    "Sample",
    "Target",
    "Cq",
    "Amp Status"
  )
  required_group_cols <- c(
    "Sample",
    "Group",
    "Group_order"
  )
  missing_qpcr <- setdiff(
    required_qpcr_cols,
    names(raw_data)
  )
  if (length(missing_qpcr) > 0) {
    cli::cli_abort(
      c(
        "Missing required columns in qPCR file.",
        "x {.val {missing_qpcr}}"
      )
    )
  }
  missing_group <- setdiff(
    required_group_cols,
    names(group_label)
  )
  if (length(missing_group) > 0) {
    cli::cli_abort(
      c(
        "Missing required columns in group label file.",
        "x {.val {missing_group}}"
      )
    )
  }
  if (!all(
    controls %in% unique(raw_data$Target)
  )) {
    missing_controls <- setdiff(
      controls,
      unique(raw_data$Target)
    )
    cli::cli_abort(
      c(
        "Control genes not found.",
        "x {.val {missing_controls}}"
      )
    )
  }
  base_group <- group_label |>
    dplyr::filter(Group_order == 1) |>
    dplyr::pull(Group) |>
    unique()
  if (length(base_group) != 1) {
    cli::cli_abort(
      "Exactly one group must have Group_order = 1."
    )
  }
  invisible(TRUE)
}

validate_palette <- function(
    palette
) {
  allowed <- c(
    "grouped",

    "khroma_bright",
    "khroma_muted",
    "khroma_vibrant",
    "khroma_light",

    "npg",
    "aaas",
    "lancet",
    "jama",
    "nejm",
    "jco",
    "okabe",
    "viridis"
  )
  if (
    is.character(palette) &&
    length(palette) == 1 &&
    (is.null(names(palette)) || !nzchar(names(palette))) &&
    !palette %in% allowed
  ) {
    stop(
      paste(
        "palette must be one of:",
        paste(
          allowed,
          collapse = ", "
        )
      )
    )
  }
}
