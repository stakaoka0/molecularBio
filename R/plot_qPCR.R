#' Analyze qPCR data
#'
#' @description
#' `plot_qPCR()` analyzes a tsv file produced by analyze_qPCR().
#'
#'
#'@param path_to_qPCR_res path to a tsv file produced by analyze_qPCR().
#'@param path_to_stats_res path to a tsv file produced by stats_qPCR().
#'@param out_dir_path path to an output directory.
#'@param controls list of internal controls.
#'      c("18S", "GAPDH", "RPL13A") (default)
#'@param plot_all,plot_internal_controls,plot_individual,save_final_table A boolean.
#'@param file_prefix prefix for files.
#'@param dpi dpi for graphs.
#'@param palette https://nanx.me/ggsci/
#'@export

plot_qPCR <- function(
    object,

    palette = "grouped",
    plot_all = TRUE,
    plot_individual = FALSE,
    plot_internal_controls = FALSE,

    show_n = TRUE,
    show_controls = TRUE,
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

    save_tiff = TRUE,
    save_svg = TRUE

) {
  stopifnot(
    inherits(object, "qpcr_analysis")
  )

  plot_data <- object$data %>%
    dplyr::mutate(
      Group = factor(
        Group,
        levels = object$group_levels
      )
    )

  group_colors <- resolve_palette(
    palette = palette,
    groups = object$group_levels,
    group_info = object$group_info
  )

  if (plot_all) {
    overview_plot <- build_overview_plot(
      object,
      colors = group_colors
    )

    save_plot(
      plot = overview_plot,
      filename_base = file.path(
        out_dir_path,
        "qPCR_plot"
      ),
      width = width,
      height = height,
      dpi = dpi
    )

  }
  if (plot_individual) {
    genes <- unique(
      object$data$Target
    )
    plots <- purrr::map(
      genes,
      ~ build_gene_plot(
        object,
        gene = .x,
        colors = group_colors,

        show_n = show_n,
        show_controls = show_controls,

        hide_ns = hide_ns,
        stats_display = stats_display,

        comparison = comparison,
        stats_group = stats_group
      )
    )
    names(plots) <- genes
    object$plots <- plots

    if (!is.null(out_dir_path)) {
      dir.create(
        out_dir_path,
        recursive = TRUE,
        showWarnings = FALSE
      )

      purrr::iwalk(
        plots,
        ~ save_plot(
          plot = .x,
          filename_base = file.path(
            out_dir_path,
            paste0(.y, "_qPCR_plot")
          ),
          width = width,
          height = height,
          dpi = dpi
        )
      )
    }

    return(
      invisible(plots)
    )
  }
}

# plot_qPCR <- function(path_to_qPCR_res,
#                       path_to_stats_res,
#                       out_dir_path,
#                       controls = c("18S", "GAPDH", "RPL13A"),
#                       plot_all = TRUE,
#                       plot_internal_controls = FALSE,
#                       plot_individual = FALSE,
#                       file_prefix = "",
#                       stats_type = "none",
#                       stats_display = "star",
#                       show_n_and_norm_genes = TRUE,
#                       stats_group = NULL,
#                       hide.ns = TRUE,
#                       dpi = 150,
#                       ylim_fold = 1.1,
#                       width = 1.5,
#                       height = 4) {
#
#   if (!stats_type %in% c("ttest", "anova", "none")) {
#     stop("stats_type must be either 'ttest' or 'anova'.")
#   }
#
#   if (!stats_display %in% c("star", "value")) {
#     stop("stats_display must be either 'star' or 'value'.")
#   }
#
#   if (!show_n_and_norm_genes %in% c(TRUE, FALSE)) {
#     stop("show_n_and_norm_genes must be either 'TRUE' or 'FALSE'.")
#   }
#
#   # Load dependent libraries----------------------------------------------------
#   suppressPackageStartupMessages(require(readr))
#   suppressPackageStartupMessages(require(dplyr))
#   suppressPackageStartupMessages(require(ggplot2))
#   suppressPackageStartupMessages(require(ggtext))
#   suppressPackageStartupMessages(require(ggpubr))
#   suppressPackageStartupMessages(require(rstatix))
#
#   # Load and modify data for plotting ------------------------------------------
#
#   # make outdir folder
#
#   dir.create(file.path(out_dir_path),
#              recursive = TRUE,
#              showWarnings = FALSE)
#
#   # Load an analyzed qPCR data.
#   plot_data <- readr::read_tsv(path_to_qPCR_res,
#                               col_types = cols())
#
#   # Define factor for ordering the group label.
#   group_fct <- plot_data %>%
#     distinct(Group, Group_order) %>%
#     arrange(Group_order) %>%
#     pull(Group)
#
#   plot_data <- plot_data %>%
#     dplyr::mutate(Group = factor(Group, level = group_fct))
#
#   color_fct <- plot_data %>%
#     dplyr::distinct(Group, Group_order, color) %>%
#     dplyr::arrange(Group) %>%
#     dplyr::pull(color) %>%
#     unique()
#
#   # Define plot theme
#   plot_theme_basic <- theme_classic() + theme(
#     title = element_text(size = 12),
#     plot.title = element_text(hjust = 0.5),
#     axis.title.x = element_text(size = 16),
#     axis.title.y = element_text(size = 16),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 12)
#   )
#
#   # Draw plots------------------------------------------------------------------
#
#   if (plot_all == TRUE) {
#
#     if (plot_internal_controls == FALSE) {
#
#       plot_data <- plot_data %>%
#         dplyr::filter(!Target %in% controls)
#
#     }
#
#     g1 <- plot_data %>%
#       dplyr::mutate(norm_exp = case_when(is.na(norm_exp) ~ -Inf,
#                                          TRUE ~ norm_exp)) %>%
#       ggplot(aes(x=Target,
#                  y=group_exp,
#                  group=factor(Group, level = group_fct))) +
#       geom_bar(position = position_dodge(.9),
#                stat = "identity",
#                color = "black", linewidth=0.5,
#                alpha = 0.8,
#                aes(fill = factor(color, level = color_fct))) +
#       geom_errorbar(aes(ymax = upper_error, ymin = lower_error),
#                     width = .2,
#                     position = position_dodge(.9)) +
#       geom_jitter(position = position_dodge(.9),
#                   fill = "black",
#                   aes(x = Target,
#                       y = norm_exp,
#                       group = interaction(Target,factor(Group, level = group_fct)))) +
#       # geom_dotplot(binaxis = "y",
#       #              stackdir = "center",
#       #              #position =  position_dodge(.9),
#       #              fill = "black",
#       #              aes(x = Target,
#       #                  y = norm_exp,
#       #                  group = interaction(Target,factor(Group, level = group_fct))),
#       #                  dotsize = 0.1) +
#       plot_theme_basic +
#       ylab(paste0("Relative Expression")) +
#       scale_fill_identity(labels = group_fct, guide = "legend") +
#       theme(legend.title = element_blank(),
#             axis.title.y = element_markdown(),
#             axis.title.x = element_blank())
#
#     if (show_n_and_norm_genes == TRUE) {
#       g1 <- g1 +
#         ylab(paste0("Relative Expression<br><span style='font-size:10pt'>(Normalized by ",
#                     paste0(controls, collapse = " "), ")</span>"))
#     }
#
#     n_targets <- plot_data %>%
#       pull(Target) %>%
#       unique() %>%
#       length()
#
#     ggsave(plot = g1,
#            file = paste0(out_dir_path, "/",
#                          file_prefix,
#                          "qPCR_plot.tiff"),
#            width = n_targets*width, height = height, dpi = dpi
#     )
#
#   }
#
#   if (plot_individual == TRUE) {
#
#     if (plot_internal_controls == FALSE) {
#
#       plot_data <- plot_data %>%
#         dplyr::filter(!Target %in% controls)
#
#     }
#
#     genes <- plot_data %>%
#       dplyr::pull(Target) %>%
#       unique()
#
#     if (stats_type %in% c("ttest", "anova")) {
#
#       stat_test <- read_tsv(path_to_stats_res, col_types =  cols())
#
#       if (!is.null(stats_group)) {
#
#         combi <- as.data.frame(stats_group) %>%
#           t() %>%
#           as_tibble(test, rownames = NULL)
#
#         colnames(combi) <- c("group1", "group2")
#
#         df1 <- stat_test %>%
#           inner_join(combi, by = c("group1", "group2"))
#
#         df2 <- stat_test %>%
#           inner_join(combi, by = c("group1" = "group2",
#                                    "group2" = "group1")
#           )
#
#         stat_test <- df1 %>%
#           bind_rows(df2)
#         # Write code to re-define y.position.
#       }
#     }
#     for (i in genes) {
#
#       max_exp <- plot_data %>%
#         dplyr::filter(Target == i) %>%
#         pull(norm_exp) %>%
#         max(na.rm = T)
#
#       g1 <- plot_data %>%
#         dplyr::filter(Target == i) %>%
#         dplyr::mutate(norm_exp = case_when(is.na(norm_exp) ~ -Inf,
#                                            TRUE ~ norm_exp)) %>%
#         ggplot(aes(x = Group,
#                    y = log2_norm_exp,
#                    group = factor(Group, level = group_fct))) +
#         geom_bar(position = "dodge",
#                  stat = "identity",
#                  color = "black", linewidth=0.5,
#                  alpha = 0.8,
#                  aes(x = Group,
#                      y = group_exp,
#                      fill = factor(color, level = color_fct))) +
#         geom_errorbar(aes(ymax = upper_error, ymin = lower_error),
#                       width=.2,
#                       position = position_dodge(.9)) +
#         geom_jitter(position = position_dodge(.9),
#                     fill = "black",
#                     aes(x = Group,
#                         y = norm_exp,
#                         group = factor(Group, level = group_fct))) +
#         # geom_dotplot(binaxis = "y",
#         #              stackdir = "center",
#         #              position = position_dodge(.9),
#         #              fill = "black",
#         #              binwidth = max_exp/30,
#         #              dotsize = 0.5,
#         #            aes(x = Group,
#         #                y = norm_exp,
#         #                group = factor(Group, level = group_fct))) +
#         coord_cartesian(ylim=c(0, max_exp*ylim_fold)) +
#         ylab(paste0("Relative Expression")) +
#         xlab(i) +
#         plot_theme_basic +
#         scale_fill_identity(labels = group_fct, guide = "legend") +
#         theme(legend.title = element_blank(),
#               axis.title.y = element_markdown(),
#               axis.text.x = element_blank())
#
#       if (stats_type == "ttest") {
#
#         stat_test_i <- stat_test %>%
#           dplyr::filter(Target == i)
#
#         if (hide.ns == TRUE &
#             TRUE %in% (0.05 > stat_test_i$p)) {
#
#           stat_test_i <- stat_test_i %>%
#             dplyr::filter(p <= 0.05)
#
#           max_exp <- stat_test_i %>%
#             pull(y.position) %>%
#             max()
#
#           if (stats_display == "value") {
#             g1 <- g1 + stat_pvalue_manual(
#               stat_test_i,
#               label = "p = {p}",
#               y.position = stat_test_i$y.position,
#               hide.ns = FALSE
#             ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
#           } else if (stats_display == "star") {
#             g1 <- g1 + stat_pvalue_manual(
#               stat_test_i,
#               label = "p.signif",
#               y.position = stat_test_i$y.position,
#               hide.ns = FALSE
#             ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
#           }
#
#         } else if (hide.ns == FALSE) {
#
#           max_exp <- stat_test_i %>%
#             pull(y.position) %>%
#             max()
#
#           if (stats_display == "value") {
#             g1 <- g1 + stat_pvalue_manual(
#               stat_test_i,
#               label = "p = {p}",
#               y.position = stat_test_i$y.position,
#               hide.ns = FALSE
#             ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
#           } else if (stats_display == "star") {
#             g1 <- g1 + stat_pvalue_manual(
#               stat_test_i,
#               label = "p.signif",
#               y.position = stat_test_i$y.position,
#               hide.ns = FALSE
#             ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
#           }
#
#
#         }
#       } else if (stats_type == "anova") {
#
#         stat_test_i <- stat_test %>%
#           dplyr::filter(Target == i)
#
#         if (hide.ns == TRUE &
#             TRUE %in% (0.05 > stat_test_i$p.adj)) {
#
#           stat_test_i <- stat_test_i %>%
#             dplyr::filter(p.adj <= 0.05)
#
#           max_exp <- stat_test_i %>%
#             pull(y.position) %>%
#             max()
#
#           if (stats_display == "value") {
#             g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
#               stat_pvalue_manual(
#                 stat_test_i, label = "{p.adj}",
#                 y.position = stat_test_i$y.position,
#                 hide.ns = FALSE
#               ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
#           } else if (stats_display == "star") {
#             g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
#               stat_pvalue_manual(
#                 stat_test_i, label = "p.adj.signif",
#                 y.position = stat_test_i$y.position,
#                 hide.ns = FALSE
#               ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
#           }
#
#
#         } else if (hide.ns == FALSE) {
#
#           max_exp <- stat_test_i %>%
#             pull(y.position) %>%
#             max()
#
#         if (stats_display == "value") {
#           g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
#             stat_pvalue_manual(
#               stat_test_i, label = "{p.adj}",
#               y.position = stat_test_i$y.position,
#               hide.ns = FALSE
#             ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
#         } else if (stats_display == "star") {
#           g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
#             stat_pvalue_manual(
#               stat_test_i, label = "p.adj.signif",
#               y.position = stat_test_i$y.position,
#               hide.ns = FALSE
#             ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
#         }
#
#
#         }
#       }
#
#       if (show_n_and_norm_genes == TRUE) {
#         g1 <- g1 +
#           ylab(paste0("Relative Expression<br><span style='font-size:10pt'>(Normalized by ",
#                       paste0(controls, collapse = " "), ")</span>")) +
#           geom_text(position = position_dodge(width = .9),
#                     aes(label = paste("n:", n), y = 0), vjust= 1.1)
#       }
#
#       n_group <- plot_data %>%
#         pull(Group) %>%
#         unique() %>%
#         length()
#       ggsave(plot = g1,
#              file = paste0(out_dir_path, "/",
#                            file_prefix,
#                            i,
#                            "_qPCR_plot.tiff"),
#              width = n_group*width, height = height, dpi = dpi
#       )
#     }
#   }
# }
