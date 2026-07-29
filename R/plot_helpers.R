#' Build individual gene plot

#'

#' @keywords internal

build_gene_plot <- function(
  object,
  gene,
  colors,
  show_n = TRUE,
  show_controls = TRUE,
  hide_ns = TRUE,
  stats_display = c("star", "value"),
  ylim_fold = 1.1,
  comparison = "all",
  stats_group = NULL

) {

  stats_display <- match.arg(stats_display)
  df <- object$data |>
    dplyr::filter(
      Target == gene
    ) |>
    dplyr::mutate(
      Group = factor(
        Group,
        levels = object$group_levels
      )
    )

  summary_df <- df |>
    dplyr::distinct(
      Group,
      group_exp,
      upper_error,
      lower_error,
      n
    )

  max_exp <- max(
    df$norm_exp,
    na.rm = TRUE
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = summary_df,
      ggplot2::aes(
        x = Group,
        y = group_exp,
        fill = Group
      ),
      color = "black",
      linewidth = 0.4,
      alpha = 0.8
    ) +
    ggplot2::geom_errorbar(
      data = summary_df,
      ggplot2::aes(
        x = Group,
        ymin = lower_error,
        ymax = upper_error
      ),
      width = 0.2
    ) +
    ggplot2::geom_jitter(
      data = df,
      ggplot2::aes(
        x = Group,
        y = norm_exp
      ),
      width = 0.1,
      height = 0,
      size = 2,
      alpha = 0.8
    ) +
    ggplot2::scale_fill_manual(
      values = colors,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      ylim = c(
        0,
        max_exp * ylim_fold
      )
    ) +
    ggplot2::labs(
      x = gene,
      y = make_y_label(
        object$controls,
        show_controls
      )
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = element_blank(),
      legend.position = "right",
    )
  if (show_n) {
    p <- p +
      ggplot2::geom_text(
        data = summary_df,
        ggplot2::aes(
          x = Group,
          y = 0,
          label = paste0("n=", n)
        ),
        vjust = 1.2,
        size = 3
      )
  }
  if (!is.null(object$stats)) {
    stat_df <- object$stats |>
      dplyr::filter(
        Target == gene
      )
    stat_df <- filter_stats_results(
      stat_df = stat_df,
      comparison = comparison,
      reference_group = object$reference_group,
      comparisons = stats_group
    )

    stat_df <- recalculate_y_position(
      stat_df,
      max_y = max(df$norm_exp, na.rm = TRUE),
      step = max(df$norm_exp, na.rm = TRUE) * 0.08
    )

    if (nrow(stat_df) > 0) {
      p_col <- NULL
      label_col <- NULL
      if ("p" %in% names(stat_df)) {
        p_col <- "p"
      }
      if ("p.adj" %in% names(stat_df)) {
        p_col <- "p.adj"
      }
      if ("signif_label" %in% names(stat_df)) {
        label_col <- "signif_label"
      } else if ("p.signif" %in% names(stat_df)) {
        label_col <- "p.signif"
      } else if ("p.adj.signif" %in% names(stat_df)) {
        label_col <- "p.adj.signif"
      }
      if (
        hide_ns &&
        !is.null(p_col)
      ) {
        stat_df <- stat_df |>
          dplyr::filter(
            .data[[p_col]] <= 0.05
          )
      }
      if (
        nrow(stat_df) > 0 &&
        !is.null(label_col)
      ) {
        if (stats_display == "value") {
          label_col <- p_col
        }
        ymax <- max(
          stat_df$y.position,
          na.rm = TRUE
        )
        p <- p +
          ggpubr::stat_pvalue_manual(
            stat_df,
            label = label_col,
            y.position = "y.position"
          ) +
          ggplot2::coord_cartesian(
            ylim = c(
              0,
              ymax * ylim_fold
            )
          )
      }
    }
  }
  p

}

build_overview_plot <- function(
  object,
  colors,
  show_controls = TRUE

) {
  summary_df <- object$data %>%
    dplyr::mutate(
      Group = factor(
        Group,
        levels = object$group_levels
      )
    ) %>%

    dplyr::distinct(
      Group,
      Target,
      group_exp,
      upper_error,
      lower_error

    )

  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = Target,
      y = group_exp,
      fill = Group
    )

  ) +

    ggplot2::geom_col(
      position = ggplot2::position_dodge(0.9),
      color = "black",
      linewidth = 0.3
    ) +

    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower_error,
        ymax = upper_error

      ),
      position = ggplot2::position_dodge(0.9),
      width = 0.2

    ) +

    ggplot2::scale_fill_manual(

      values = colors,

      drop = FALSE

    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = NULL,
      y = make_y_label(
        object$controls,
        show_controls
      ),
      fill = NULL
    )

}

recalculate_y_position <- function(
    stat_df,
    max_y,
    step = 0.1
) {
  if (nrow(stat_df) == 0) {
    return(stat_df)
  }

  stat_df <- stat_df |>
    dplyr::arrange(group1, group2)
  stat_df$y.position <-
    max_y + step * seq_len(nrow(stat_df))
  stat_df
}
