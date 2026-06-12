#' Wrapper function for qPCR analysis and plotting graphs
#'
#' @description
#' `qPCR()` Wrapper function for qPCR analysis and plotting graphs
#'
#'@param path_to_qPCR_res path to a csv file produced by QuantStudio.
#'@param path_to_group_label path to a csv file including meta data.
#'@param out_dir_path path to an output directory.
#'@param controls list of internal controls.
#'      c("18S", "GAPDH", "RPL13A") (default)
#'@param save_final_table A boolean.
#'@param file_prefix prefix for files.
#'@export
#'
#'
#'

qPCR <- function(path_to_qPCR_res,
                 path_to_group_label,
                 out_dir_path,
                 controls = c("18S", "GAPDH", "RPL13A"),
                 save_final_table = TRUE,
                 file_prefix = "",
                 stats_type = "none",
                 stats_value = "log2_norm_exp",
                 stats_display = "star",
                 show_n_and_norm_genes = TRUE,
                 paired = F,
                 var.equal = TRUE,
                 alternative = "two.sided",
                 stats_group = NULL,
                 plot_all = TRUE,
                 plot_internal_controls = FALSE,
                 plot_individual = FALSE,
                 hide.ns = TRUE,
                 dpi = 150,
                 ylim_fold = 1.1,
                 width = 1.5,
                 height = 4) {

  analyze_qPCR(path_to_qPCR_res = path_to_qPCR_res,
               path_to_group_label = path_to_group_label,
               out_dir_path = out_dir_path,
               controls = controls,
               save_final_table = save_final_table,
               file_prefix = file_prefix)


  if (stats_type != "none") {
    stats_qPCR(path_to_qPCR_res = paste0(out_dir_path,
                                         "/",
                                         file_prefix,
                                         "analyzed_qPCR_res.tsv"),
               out_dir_path = out_dir_path,
               controls = controls,
               stats_type = stats_type,
               stats_value = stats_value,
               paired = paired,
               var.equal = var.equal,
               alternative = alternative,
               stats_group = stats_group,
               file_prefix = file_prefix)

    if (stats_type == "ttest") {
      path_to_stats_res <- paste0(out_dir_path,
                                 "/",
                                 file_prefix,
                                 "qPCR_ttest_result.tsv")
    } else if (stats_type == "anova") {
      path_to_stats_res <- paste0(out_dir_path,
                                 "/",
                                 file_prefix,
                                 "qPCR_anova_tukey_result.tsv")
    }

  }

  plot_qPCR(path_to_qPCR_res = paste0(out_dir_path,
                                      "/",
                                      file_prefix,
                                      "analyzed_qPCR_res.tsv"),
            path_to_stats_res = path_to_stats_res,
            out_dir_path = out_dir_path,
            controls = controls,
            plot_all = plot_all,
            plot_internal_controls = plot_internal_controls,
            plot_individual = plot_individual,
            file_prefix = file_prefix,
            stats_type = stats_type,
            stats_display = stats_display,
            stats_group = stats_group,
            show_n_and_norm_genes = show_n_and_norm_genes,
            hide.ns = hide.ns,
            dpi = dpi,
            ylim_fold = ylim_fold,
            width = width,
            height = height)

}



