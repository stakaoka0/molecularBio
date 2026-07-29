#' Plot qPCR results
#'
#' Creates an overview plot, individual target plots, or both. Statistical
#' annotations are added when `object` was processed by [stats_qPCR()].
#'
#' @param object A `qpcr_analysis` or `qpcr_stats` object.
#' @param palette A supported palette name or a character vector of colors.
#' @param plot_all Logical; create the overview plot.
#' @param plot_individual Logical; create one plot per target.
#' @param plot_internal_controls Logical; include reference genes.
#' @param show_n Logical; show sample size on individual plots.
#' @param show_controls Logical; list reference genes in the y-axis label.
#' @param jitter Logical; horizontally jitter sample points on individual plots.
#' @param jitter_width Non-negative numeric jitter width.
#' @param jitter_seed Integer seed used to make jitter positions reproducible.
#' @param comparison Statistical comparisons to show.
#' @param stats_group List of two-element character vectors used when
#'   `comparison = "custom"`.
#' @param hide_ns Logical; hide nonsignificant comparisons.
#' @param stats_display Display significance as `"star"` or `"value"`.
#' @param out_dir_path Directory for saved plots.
#' @param width,height Plot dimensions in inches.
#' @param dpi Raster resolution.
#' @param save_tiff,save_svg Logical; save plots in the corresponding format.
#'
#' @return Invisibly, a named list containing the generated plots.
#' @export

plot_qPCR <- function(
    object,

    palette = "grouped",
    plot_all = TRUE,
    plot_individual = FALSE,
    plot_internal_controls = FALSE,

    show_n = TRUE,
    show_controls = TRUE,
    jitter = TRUE,
    jitter_width = 0.1,
    jitter_seed = 1L,
    comparison = c(
      "all",
      "control",
      "custom",
      "significant"
    ),
    stats_group = NULL,
    hide_ns = FALSE,
    stats_display = c("star", "value"),

    out_dir_path = NULL,
    width = 1.5,
    height = 4,
    dpi = 300,

    save_tiff = FALSE,
    save_svg = FALSE

) {
  if (!inherits(object, "qpcr_analysis")) {
    cli::cli_abort("{.arg object} must be a `qpcr_analysis` object.")
  }
  if (!is.logical(jitter) || length(jitter) != 1L || is.na(jitter)) {
    cli::cli_abort("{.arg jitter} must be `TRUE` or `FALSE`.")
  }
  if (!is.numeric(jitter_width) || length(jitter_width) != 1L ||
      is.na(jitter_width) || jitter_width < 0) {
    cli::cli_abort(
      "{.arg jitter_width} must be one non-negative number."
    )
  }
  if (!is.numeric(jitter_seed) || length(jitter_seed) != 1L ||
      is.na(jitter_seed) || jitter_seed %% 1 != 0) {
    cli::cli_abort("{.arg jitter_seed} must be one integer.")
  }
  jitter_seed <- as.integer(jitter_seed)
  comparison <- match.arg(comparison)
  stats_display <- match.arg(stats_display)
  validate_palette(palette)
  if (!is.logical(save_tiff) || length(save_tiff) != 1L ||
      is.na(save_tiff) ||
      !is.logical(save_svg) || length(save_svg) != 1L ||
      is.na(save_svg)) {
    cli::cli_abort("{.arg save_tiff} and {.arg save_svg} must be logical.")
  }
  save_requested <- save_tiff || save_svg
  if (save_requested && is.null(out_dir_path)) {
    cli::cli_abort(
      "{.arg out_dir_path} is required when plot saving is enabled."
    )
  }
  if (comparison == "custom" && is.null(stats_group)) {
    cli::cli_abort(
      "{.arg stats_group} is required when {.arg comparison} is `\"custom\"`."
    )
  }

  group_colors <- resolve_palette(
    palette = palette,
    groups = object$group_levels,
    group_info = object$group_info
  )
  plots <- list()

  if (plot_all) {
    overview_plot <- build_overview_plot(
      object,
      colors = group_colors,
      show_controls = show_controls,
      include_internal_controls = plot_internal_controls
    )
    plots$overview <- overview_plot

    if (save_requested) {
      dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)
      save_plot(
        plot = overview_plot,
        filename_base = file.path(out_dir_path, "qPCR_plot"),
        width = width,
        height = height,
        dpi = dpi,
        save_tiff = save_tiff,
        save_svg = save_svg
      )
    }
  }
  if (plot_individual) {
    genes <- unique(
      object$data$Target[
        plot_internal_controls |
          !object$data$Target %in% object$controls
      ]
    )
    individual_plots <- purrr::map(
      genes,
      ~ build_gene_plot(
        object,
        gene = .x,
        colors = group_colors,

        show_n = show_n,
        show_controls = show_controls,
        jitter = jitter,
        jitter_width = jitter_width,
        jitter_seed = jitter_seed,

        hide_ns = hide_ns,
        stats_display = stats_display,

        comparison = comparison,
        stats_group = stats_group
      )
    )
    names(individual_plots) <- genes
    plots$individual <- individual_plots

    if (save_requested) {
      purrr::iwalk(
        individual_plots,
        ~ save_plot(
          plot = .x,
          filename_base = file.path(
            out_dir_path,
            paste0(.y, "_qPCR_plot")
          ),
          width = width,
          height = height,
          dpi = dpi,
          save_tiff = save_tiff,
          save_svg = save_svg
        )
      )
    }
  }
  invisible(plots)
}
