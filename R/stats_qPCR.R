#' Perform statistical analysis on qPCR data
#'
#' Runs gene-wise pairwise t-tests or one-way ANOVA followed by Tukey tests.
#' Reference genes are excluded from testing.
#'
#' @param object A `qpcr_analysis` object returned by [analyze_qPCR()].
#' @param stats_type Statistical method: `"ttest"` or `"anova"`.
#' @param stats_value Expression column to test.
#' @param paired Logical; use paired t-tests.
#' @param var.equal Logical; assume equal variances for t-tests.
#' @param alternative Alternative hypothesis for t-tests.
#' @param out_dir_path Output directory used when `save_results = TRUE`.
#' @param save_results Logical; save results as `qPCR_stats.tsv`.
#'
#' @return The input object with a `stats` component and class `qpcr_stats`.
#' @export

stats_qPCR <- function(
  object,
  stats_type = c("ttest", "anova"),
  stats_value = c("log2_norm_exp", "norm_exp"),
  paired = FALSE,
  var.equal = TRUE,
  alternative = c("two.sided", "less", "greater"),
  out_dir_path = NULL,
  save_results = FALSE
) {
  if (!inherits(object, "qpcr_analysis")) {
    cli::cli_abort("{.arg object} must be a `qpcr_analysis` object.")
  }
  stats_type <- match.arg(stats_type)
  stats_value <- match.arg(stats_value)
  alternative <- match.arg(alternative)
  if (!is.logical(save_results) || length(save_results) != 1L ||
      is.na(save_results)) {
    cli::cli_abort("{.arg save_results} must be `TRUE` or `FALSE`.")
  }
  if (save_results && is.null(out_dir_path)) {
    cli::cli_abort(
      "{.arg out_dir_path} is required when {.arg save_results} is `TRUE`."
    )
  }
  formula <- stats::as.formula(
    paste(stats_value, "~ Group")
  )
  results <- object$data |>
    dplyr::filter(
      !Target %in% object$controls
    ) |>
    dplyr::group_split(Target) |>
    purrr::map_dfr(
      function(df) {
        gene <- unique(df$Target)
        tryCatch(
          run_gene_statistics(
            df,
            formula = formula,
            stats_type = stats_type,
            paired = paired,
            var.equal = var.equal,
            alternative = alternative
          ),
          error = function(error) {
            cli::cli_warn(
              "Skipping {gene}: {conditionMessage(error)}"
            )
            tibble::tibble()
          }
        )
      }
    )
  object$stats <- results

  if (save_results) {
    dir.create(
      out_dir_path,
      recursive = TRUE,
      showWarnings = FALSE
    )
    readr::write_tsv(
      results,
      file.path(
        out_dir_path,
        "qPCR_stats.tsv"
      )
    )
  }

  class(object) <- c(
    "qpcr_stats",
    class(object)
  )

  object
}
