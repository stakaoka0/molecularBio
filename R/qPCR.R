#' Run the complete qPCR workflow
#'
#' Runs [analyze_qPCR()], optionally [stats_qPCR()], and then [plot_qPCR()].
#'
#' @inheritParams analyze_qPCR
#' @inheritParams stats_qPCR
#' @inheritParams plot_qPCR
#' @param stats_type Statistical method, or `"none"` to skip statistics.
#' @param save_analysis Logical; save the analyzed expression table.
#' @param save_stats Logical; save the statistical results.
#' @param save_plots Logical; save both TIFF and SVG plots.
#'
#' @return Invisibly, a `qpcr_analysis` or `qpcr_stats` object.
#' @export
qPCR <- function(

  path_to_qPCR_res,
  path_to_group_label,
  out_dir_path,
  controls = c("18S", "GAPDH", "RPL13A"),
  technical_replicates = c("average", "error"),
  technical_sd_threshold = 0.2,
  stats_type = c("none", "ttest", "anova"),
  stats_value = c("log2_norm_exp", "norm_exp"),
  paired = FALSE,
  var.equal = TRUE,
  alternative = c("two.sided", "less", "greater"),
  comparison = c(
    "all",
    "control",
    "custom",
    "significant"
  ),
  stats_group = NULL,
  plot_all = TRUE,
  plot_individual = FALSE,
  plot_internal_controls = FALSE,
  hide_ns = TRUE,
  width = 1.5,
  height = 4,
  dpi = 300,
  save_analysis = TRUE,
  save_stats = TRUE,
  save_plots = TRUE,
  palette = "grouped",
  stats_display = c("star", "value"),
  show_n = TRUE,
  show_controls = TRUE,
  jitter = TRUE,
  jitter_width = 0.1,
  jitter_seed = 1L

) {
  stats_type <- match.arg(stats_type)
  stats_value <- match.arg(stats_value)
  alternative <- match.arg(alternative)
  comparison <- match.arg(comparison)
  stats_display <- match.arg(stats_display)

  res <- analyze_qPCR(

    path_to_qPCR_res = path_to_qPCR_res,
    path_to_group_label = path_to_group_label,
    controls = controls,
    technical_replicates = technical_replicates,
    technical_sd_threshold = technical_sd_threshold,
    out_dir_path = out_dir_path,
    save_results = save_analysis

  )

  if (stats_type != "none") {

    res <- stats_qPCR(
      res,
      stats_type = stats_type,
      stats_value = stats_value,
      paired = paired,
      var.equal = var.equal,
      alternative = alternative,
      out_dir_path = out_dir_path,
      save_results = save_stats

    )

  }

  plot_qPCR(

    object = res,
    palette = palette,

    plot_all = plot_all,
    plot_individual = plot_individual,
    plot_internal_controls = plot_internal_controls,
    show_n = show_n,
    show_controls = show_controls,
    jitter = jitter,
    jitter_width = jitter_width,
    jitter_seed = jitter_seed,

    comparison = comparison,
    stats_group = stats_group,
    hide_ns = hide_ns,
    stats_display = stats_display,

    out_dir_path = out_dir_path,

    width = width,
    height = height,
    dpi = dpi,
    save_tiff = save_plots,
    save_svg = save_plots

  )

  invisible(res)

}
