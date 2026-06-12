#' Analyze qPCR data
#'
#' @description
#' `plot_qPCR_treatment()` analyzes a tsv file produced by analyze_qPCR_treatment().
#'
#'
#'@param path_to_qPCR_res path to a tsv file produced by analyze_qPCR_treatment().
#'@param path_to_stats_res path to a tsv file produced by stats_qPCR_treatment().
#'@param out_dir_path path to an output directory.
#'@param controls list of internal controls.
#'      c("18S", "GAPDH", "RPL13A") (default)
#'@param plot_all,plot_internal_controls,plot_individual,save_final_table A boolean.
#'@param file_prefix prefix for files.
#'@param dpi dpi for graphs.
#'@export

plot_qPCR_treatment <- function(path_to_qPCR_res,
                      path_to_stats_res,
                      out_dir_path,
                      controls = c("18S", "GAPDH", "RPL13A"),
                      plot_all = TRUE,
                      plot_internal_controls = FALSE,
                      plot_individual = FALSE,
                      file_prefix = "",
                      stats_type = "none",
                      show_norm_genes = TRUE,
                      stats_group = NULL,
                      hide.ns = TRUE,
                      dpi = 150,
                      ylim_fold = 1.1,
                      width = 2.5,
                      height = 4) {

  if (!stats_type %in% c("ttest", "anova", "none")) {
    stop("stats_type must be either 'ttest' or 'anova'.")
  }

  if (!show_norm_genes %in% c(TRUE, FALSE)) {
    stop("show_norm_genes must be either 'TRUE' or 'FALSE'.")
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

  # Load an analyzed qPCR data.
  plot_data <- readr::read_tsv(path_to_qPCR_res,
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

  conc_list <- plot_data %>%
    pull(treatment_conc) %>%
    sort() %>%
    unique()

  conc_list <- factor(conc_list, levels = conc_list)

  plot_data <- plot_data %>%
    dplyr::mutate(treatment_conc = factor(treatment_conc, levels=conc_list))

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

  if (plot_all == TRUE) {

    if (plot_internal_controls == FALSE) {

      plot_data <- plot_data %>%
        dplyr::filter(!Target %in% controls)

    }

    g1 <- plot_data %>%
      dplyr::mutate(norm_exp = case_when(is.na(norm_exp) ~ -Inf,
                                         TRUE ~ norm_exp)) %>%
      ggplot(aes(x = treatment_conc,
                 y = group_exp)) +
      facet_wrap(~Target, scales = "free") +
      geom_line(aes(color = factor(color, level = color_fct))) +
      expand_limits(y=0) +
      geom_errorbar(aes(ymax = upper_error, ymin = lower_error),
                    width = .2,
                    position = position_dodge(.9)) +
      geom_dotplot(binaxis = "y",
                   stackdir = "center",
                   #position =  position_dodge(.9),
                   aes(x = treatment_conc,
                       y = norm_exp,
                       group = interaction(Target,factor(treatment_conc))),
                       dotsize = 0.5) +
      geom_hline(yintercept = 0.5, linetype="dotted") +
      ylab(paste0("Relative Expression")) +
      plot_theme_basic +
      scale_color_identity(labels = group_fct, guide = "legend") +
      scale_fill_identity(labels = group_fct, guide = "legend") +
      scale_x_discrete(breaks = conc_list) +
      theme(legend.title = element_blank(),
            axis.title.y = element_markdown(),
            axis.title.x = element_blank())

    if (show_norm_genes == TRUE) {
      g1 <- g1 +
        ylab(paste0("Relative Expression<br><span style='font-size:10pt'>(Normalized by ",
                    paste0(controls, collapse = " "), ")</span>"))
    }

    n_targets <- plot_data %>%
      pull(Target) %>%
      unique() %>%
      length()

    ggsave(plot = g1,
           file = paste0(out_dir_path, "/",
                         file_prefix,
                         "qPCR_plot.tiff"),
           width = n_targets*width, height = height, dpi = dpi
    )

  }

  if (plot_individual == TRUE) {

    if (plot_internal_controls == FALSE) {

      plot_data <- plot_data %>%
        dplyr::filter(!Target %in% controls)

    }

    genes <- plot_data %>%
      dplyr::pull(Target) %>%
      unique()

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
    for (i in genes) {

      max_exp <- plot_data %>%
        dplyr::filter(Target == i) %>%
        pull(norm_exp) %>%
        max(na.rm = T)

      g1 <- plot_data %>%
        dplyr::filter(Target == i) %>%
        dplyr::mutate(norm_exp = case_when(is.na(norm_exp) ~ -Inf,
                                           TRUE ~ norm_exp)) %>%
        ggplot(aes(x = treatment_conc,
                   y = group_exp,
                   shape = Group)) +
        geom_line(aes(group = Group, color = factor(color, level = color_fct))) +
        # geom_errorbar(aes(ymax = upper_error,
        #                   ymin = lower_error,
        #                   color = factor(color, level = color_fct)),
        #               width = .2) +
        # geom_text(position = position_dodge(width = 1.5),
        #           aes(label = paste0("n=", n),
        #               y = 0,
        #               color = factor(color, level = color_fct)),
        #           vjust= 1.1) +
        geom_dotplot(binaxis = "y",
                     stackdir = "center",
                     #position =  position_dodge(.2),
                     fill = "white",
                     aes(x = treatment_conc,
                         y = norm_exp,
                         group = interaction(Group, factor(treatment_conc)),
                         color = factor(color, level = color_fct)),
                     dotsize = 0.5) +
        coord_cartesian(ylim=c(0, max_exp*ylim_fold)) +
        ylab(paste0("Relative Expression")) +
        ggtitle(i) + xlab("Concentration of GW7647 (µM)") +
        plot_theme_basic +
        scale_color_identity(labels = group_fct, guide = "legend") +
        scale_fill_identity(labels = group_fct, guide = "legend") +
        scale_x_discrete(breaks = conc_list) +
        theme(legend.title = element_blank(),
              axis.title.y = element_markdown())

      if (stats_type == "ttest") {

        stat_test_i <- stat_test %>%
          dplyr::filter(Target == i) %>%
          dplyr::mutate(Group = group_fct[1])
          #group1 = 0,
        #                 #group2 = conc_list,
        #                 )

        if (hide.ns == TRUE &
            TRUE %in% (0.05 > stat_test_i$p)) {

          stat_test_i <- stat_test_i %>%
            dplyr::filter(p <= 0.05)

          max_exp <- stat_test_i %>%
            pull(y.position) %>%
            max()

          g1 <- g1 + stat_pvalue_manual(
            stat_test_i,
            label = "p = {p}",
            y.position = stat_test_i$y.position,
            hide.ns = FALSE,
            remove.bracket = TRUE
          ) + coord_cartesian(ylim = c(0, max_exp*ylim_fold))

        } else if (hide.ns == FALSE) {

          max_exp <- stat_test_i %>%
            pull(y.position) %>%
            max()

          g1 <- g1 + stat_pvalue_manual(
            stat_test_i,
            label = "p = {p}",
            y.position = stat_test_i$y.position * ylim_fold,
            hide.ns = FALSE,
            remove.bracket = TRUE
          ) + coord_cartesian(ylim = c(0, max_exp*ylim_fold))
        }
      } else if (stats_type == "anova") {

        cat("Anova for treatment chase experiments is not implemented.\n")

        # stat_test_i <- stat_test %>%
        #   dplyr::filter(Target == i)
        #
        # if (hide.ns == TRUE &
        #     TRUE %in% (0.05 > stat_test_i$p.adj)) {
        #
        #   stat_test_i <- stat_test_i %>%
        #     dplyr::filter(p.adj <= 0.05)
        #
        #   max_exp <- stat_test_i %>%
        #     pull(y.position) %>%
        #     max()
        #
        #   g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
        #     stat_pvalue_manual(
        #       stat_test_i, label = "{p.adj}",
        #       y.position = stat_test_i$y.position,
        #       hide.ns = FALSE
        #     ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
        #
        # } else if (hide.ns == FALSE) {
        #
        #   max_exp <- stat_test_i %>%
        #     pull(y.position) %>%
        #     max()
        #
        # g1 <- g1 + stat_compare_means(method = "anova", vjust = -3) +
        #   stat_pvalue_manual(
        #     stat_test_i, label = "{p.adj}",
        #     y.position = stat_test_i$y.position,
        #     hide.ns = FALSE
        #   ) + coord_cartesian(ylim=c(0, max_exp*ylim_fold))
        #
        # }
      }

      if (show_norm_genes == TRUE) {
        g1 <- g1 +
          ylab(paste0("Relative Expression<br><span style='font-size:10pt'>(Normalized by ",
                      paste0(controls, collapse = " "), ")</span>"))
      }

      n_group <- plot_data %>%
        pull(Group) %>%
        unique() %>%
        length()

      ggsave(plot = g1,
             file = paste0(out_dir_path, "/",
                           file_prefix,
                           i,
                           "_qPCR_plot.tiff"),
             width = n_group*width, height = height, dpi = dpi
      )
    }
  }
}
