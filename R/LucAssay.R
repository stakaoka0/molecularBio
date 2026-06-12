#' Wrapper function for qPCR analysis and plotting graphs
#'
#' @description
#' `LucAssay()` Wrapper function for Luciferase Assay analysis and plotting graphs
#'
#'@param path_to_firefly_res path to a csv file produced by QuantStudio.
#'@param path_to_renilla_res path to a csv file including meta data.
#'@param out_dir_path path to an output directory.
#'@param path_to_group_label
#'@param save_final_table A boolean.
#'@param file_prefix prefix for files.
#'@export
#'
#'
#'

LucAssay <- function(path_to_firefly_res,
                 path_to_renilla_res,
                 path_to_group_label,
                 out_dir_path,
                 save_final_table = TRUE,
                 file_prefix = "",
                 stats_type = "none",
                 stats_display = "star",
                 show_n= TRUE,
                 paired = F,
                 var.equal = TRUE,
                 alternative = "two.sided",
                 stats_group = NULL,
                 hide.ns = TRUE,
                 dpi = 150,
                 ylim_fold = 1.1,
                 width = 1.5,
                 height = 4) {

  analyze_LucAssay(path_to_firefly_res = path_to_firefly_res,
                   path_to_renilla_res = path_to_renilla_res,
                   path_to_group_label = path_to_group_label,
                   out_dir_path = out_dir_path,
                   save_final_table = save_final_table,
                   file_prefix = file_prefix)


  if (stats_type != "none") {
    stats_LucAssay(path_to_LucAssay_res = paste0(out_dir_path,
                                         "/",
                                         file_prefix,
                                         "analyzed_LucAssay_res.tsv"),
                   out_dir_path = out_dir_path,
                   stats_type = stats_type,
                   paired = paired,
                   var.equal = var.equal,
                   alternative = alternative,
                   stats_group = stats_group,
                   file_prefix = file_prefix)

    if (stats_type == "ttest") {
      path_to_stats_res <- paste0(out_dir_path,
                                 "/",
                                 file_prefix,
                                 "LucAssay_ttest_result.tsv")
    } else if (stats_type == "anova") {
      path_to_stats_res <- paste0(out_dir_path,
                                 "/",
                                 file_prefix,
                                 "LucAssay_anova_tukey_result.tsv")
    }

  }

  plot_LucAssay(path_to_LucAssay_res = paste0(out_dir_path,
                                      "/",
                                      file_prefix,
                                      "analyzed_LucAssay_res.tsv"),
            path_to_stats_res = path_to_stats_res,
            out_dir_path = out_dir_path,
            file_prefix = file_prefix,
            stats_type = stats_type,
            stats_display = stats_display,
            stats_group = stats_group,
            show_n = show_n,
            hide.ns = hide.ns,
            dpi = dpi,
            ylim_fold = ylim_fold,
            width = width,
            height = height)

}



