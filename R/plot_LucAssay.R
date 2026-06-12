#' Analyze LucAssay data
#'
#' @description
#' `plot_LucAssay()` analyzes a tsv file produced by analyze_LucAssay().
#'
#'
#'@param path_to_LucAssay_res path to a tsv file produced by analyze_LucAssay().
#'@param path_to_stats_res path to a tsv file produced by stats_LucAssay().
#'@param out_dir_path path to an output directory.
#'@param plot_all,plot_internal_controls,plot_individual,save_final_table A boolean.
#'@param file_prefix prefix for files.
#'@param dpi dpi for graphs.
#'@export

plot_LucAssay <- function(path_to_LucAssay_res,
                      path_to_stats_res,
                      out_dir_path,
                      show_n = TRUE,
                      file_prefix = "",
                      stats_type = "none",
                      stats_display = "star",
                      stats_group = NULL,
                      hide.ns = TRUE,
                      dpi = 150,
                      ylim_fold = 1.1,
                      width = 1.5,
                      height = 4) {

  if (!stats_type %in% c("ttest", "anova", "none")) {
    stop("stats_type must be either 'ttest' or 'anova'.")
  }

  if (!stats_display %in% c("star", "value")) {
    stop("stats_display must be either 'star' or 'value'.")
  }

  if (!show_n %in% c(TRUE, FALSE)) {
    stop("show_n must be either 'TRUE' or 'FALSE'.")
  }

  # Load dependent libraries----------------------------------------------------
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(ggtext))
  suppressPackageStartupMessages(require(ggpubr))
  suppressPackageStartupMessages(require(rstatix))

  # Load and modify data for plotting ------------------------------------------

  # make outdir folder

  dir.create(file.path(out_dir_path),
             recursive = TRUE,
             showWarnings = FALSE)

  # Load an analyzed LucAssay data.
  plot_data <- readr::read_tsv(path_to_LucAssay_res,
                              col_types = cols())

  # Define factor for ordering the group label.
  group_fct <- plot_data %>%
    distinct(Group, Group_order) %>%
    arrange(Group_order) %>%
    pull(Group)

  plot_data <- plot_data %>%
    dplyr::mutate(Group = factor(Group, level = group_fct))

  color_fct <- plot_data %>%
    dplyr::distinct(Group, Group_order, color) %>%
    dplyr::arrange(Group) %>%
    dplyr::pull(color) %>%
    unique()

  # Define plot theme
  plot_theme_basic <- theme_classic() + theme(
    title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

  # Draw plots------------------------------------------------------------------


  for (i in c("Firefly_Luc_Signal", "Renilla_Luc_Signal", "FR_Ratio")) {

    if (stats_type %in% c("ttest", "anova")) {

      stat_test <- read_tsv(path_to_stats_res, col_types =  cols())

      if (!is.null(stats_group)) {

        combi <- as.data.frame(stats_group) %>%
          t() %>%
          as_tibble(test, rownames = NULL)

        colnames(combi) <- c("group1", "group2")

        df1 <- stat_test %>%
          inner_join(combi, by = c("group1", "group2"))

        df2 <- stat_test %>%
          inner_join(combi, by = c("group1" = "group2",
                                   "group2" = "group1")
          )

        stat_test <- df1 %>%
          bind_rows(df2)
        # Write code to re-define y.position.
      }
    }

      max_sig <- plot_data %>%
        pull(!!i) %>%
        max(na.rm = T)

      g1 <- plot_data %>%
        ggplot(aes(x = Group,
                   y = eval(parse(text = i)),
                   group = factor(Group, level = group_fct))) +
        geom_bar(position = position_dodge(.9),
                 stat = "identity",
                 color = "black", linewidth=0.5,
                 alpha = 0.8,
                 aes(x = Group,
                     y = eval(parse(text = paste0("Ave_", i))),
                     fill = factor(color, level = color_fct))) +
        geom_dotplot(binaxis = "y",
                     stackdir = "center",
                     position =  position_dodge(.9),
                     fill = "black",
                     aes(x = Group,
                         y = eval(parse(text = i)),
                         group = factor(Group, level = group_fct)),
                     dotsize = 0.3) +

        coord_cartesian(ylim=c(0, max_sig*ylim_fold)) +
        xlab(i) +
        plot_theme_basic +
        scale_fill_identity(labels = group_fct, guide = "legend") +
        theme(legend.title = element_blank(),
              axis.title.y = element_markdown(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank())

      if (i == "FR_Ratio") {
        g1 <- g1 + geom_errorbar(aes(ymax = Ave_FR_Ratio+SEM_FR_Ratio,
                                     ymin = Ave_FR_Ratio-SEM_FR_Ratio),
                                 width = .2,
                                 position = position_dodge(.9)) +
          ylab(paste0("Ratio (Firefly RLU/Renilla RLU)"))
      } else if (i == "Firefly_Luc_Signal") {

        g1 <- g1 + geom_errorbar(aes(ymax = Ave_Firefly_Luc_Signal+SEM_FL,
                                     ymin = Ave_Firefly_Luc_Signal-SEM_FL),
                                 width = .2,
                                 position = position_dodge(.9)) +
          ylab(paste0("Firefly Luciferase RLU"))
      } else if (i == "Renilla_Luc_Signal") {
        g1 <- g1 + geom_errorbar(aes(ymax = Ave_Renilla_Luc_Signal+SEM_RL,
                                     ymin = Ave_Renilla_Luc_Signal-SEM_RL),
                                 width = .2,
                                 position = position_dodge(.9)) +
          ylab(paste0("Renilla Luciferase RLU"))
      }

      if (stats_type == "ttest") {

        stat_test_i <- stat_test %>%
          dplyr::filter(variables == i)

        if (hide.ns == TRUE &
            TRUE %in% (0.05 > stat_test_i$p)) {

          stat_test_i <- stat_test_i %>%
            dplyr::filter(p <= 0.05)

          max_sig <- stat_test_i %>%
            pull(y.position) %>%
            max()

          if (stats_display == "value") {
            g1 <- g1 + stat_pvalue_manual(
              stat_test_i,
              label = "p = {p}",
              y.position = stat_test_i$y.position,
              hide.ns = FALSE
            ) + coord_cartesian(ylim=c(0, max_sig*ylim_fold))
          } else if (stats_display == "star") {
            g1 <- g1 + stat_pvalue_manual(
              stat_test_i,
              label = "p.signif",
              y.position = stat_test_i$y.position,
              hide.ns = FALSE
            ) + coord_cartesian(ylim=c(0, max_sig*ylim_fold))
          }

        } else if (hide.ns == FALSE) {

          max_sig <- stat_test_i %>%
            pull(y.position) %>%
            max()

          if (stats_display == "value") {
            g1 <- g1 + stat_pvalue_manual(
              stat_test_i,
              label = "p = {p}",
              y.position = stat_test_i$y.position,
              hide.ns = FALSE
            ) + coord_cartesian(ylim=c(0, max_sig*ylim_fold))
          } else if (stats_display == "star") {
            g1 <- g1 + stat_pvalue_manual(
              stat_test_i,
              label = "p.signif",
              y.position = stat_test_i$y.position,
              hide.ns = FALSE
            ) + coord_cartesian(ylim=c(0, max_sig*ylim_fold))
          }


        }
      } else if (stats_type == "anova") {

        stat_test_i <- stat_test %>%
          dplyr::filter(variables == i)

        if (hide.ns == TRUE &
            TRUE %in% (0.05 > stat_test_i$p.adj)) {

          stat_test_i <- stat_test_i %>%
            dplyr::filter(p.adj <= 0.05)

          max_sig <- stat_test_i %>%
            pull(y.position) %>%
            max()

          if (stats_display == "value") {
            g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
              stat_pvalue_manual(
                stat_test_i, label = "{p.adj}",
                y.position = stat_test_i$y.position,
                hide.ns = FALSE
              ) + coord_cartesian(ylim=c(0, max_sig*ylim_fold))
          } else if (stats_display == "star") {
            g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
              stat_pvalue_manual(
                stat_test_i, label = "p.adj.signif",
                y.position = stat_test_i$y.position,
                hide.ns = FALSE
              ) + coord_cartesian(ylim=c(0, max_sig*ylim_fold))
          }


        } else if (hide.ns == FALSE) {

          max_sig <- stat_test_i %>%
            pull(y.position) %>%
            max()

        if (stats_display == "value") {
          g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
            stat_pvalue_manual(
              stat_test_i, label = "{p.adj}",
              y.position = stat_test_i$y.position,
              hide.ns = FALSE
            ) + coord_cartesian(ylim=c(0, max_sig*ylim_fold))
        } else if (stats_display == "star") {
          g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
            stat_pvalue_manual(
              stat_test_i, label = "p.adj.signif",
              y.position = stat_test_i$y.position,
              hide.ns = FALSE
            ) + coord_cartesian(ylim=c(0, max_sig*ylim_fold))
        }


        }
      }

      if (show_n == TRUE) {
        g1 <- g1 + geom_text(position = position_dodge(width = .9),
                             aes(label = paste("n:", n), y = 0), vjust= 1.1)
      }

      n_group <- plot_data %>%
        pull(Group) %>%
        unique() %>%
        length()
      ggsave(plot = g1,
             file = paste0(out_dir_path, "/",
                           file_prefix,
                           i,
                           "_LucAssay_plot.tiff"),
             width = n_group*width, height = height, dpi = dpi
      )
    }

}
