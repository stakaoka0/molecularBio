#' Perform statistical analysis on qPCR data
#'
#' @description
#' `stats_qPCR()` analyzes a tsv file produced by analyze_qPCR(). \n
#' This function implements the method described in this reference. \n
#' Reference https://www.cell.com/trends/biotechnology/fulltext/S0167-7799(18)30342-1
#'
#'@param path_to_qPCR_res path to a tsv file produced by analyze_qPCR().
#'@param out_dir_path path to an output directory.
#'@param controls list of internal controls.
#'      c("18S", "GAPDH", "RPL13A") (default)
#'@param file_prefix prefix for files.
#'@export

stats_qPCR <- function(
  object,
  stats_type = c("ttest", "anova"),
  stats_value = "log2_norm_exp",
  stats_group = NULL,
  paired = FALSE,
  var.equal = TRUE,
  alternative = "two.sided",
  out_dir_path = NULL,
  save_results = TRUE
) {
  stopifnot(
    inherits(object, "qpcr_analysis")
  )
  stats_type <- match.arg(
    stats_type
  )
  formula <- stats::as.formula(
    paste(stats_value, "~ Group")
  )
  results <- object$data |>
    dplyr::filter(
      !Target %in% object$controls
    ) |>
    dplyr::group_split(Target) |>
    purrr::map_dfr(
      run_gene_statistics,
      formula = formula,
      stats_type = stats_type,
      paired = paired,
      var.equal = var.equal,
      alternative = alternative
    )
  object$stats <- results

  if (
    !is.null(out_dir_path)
  ) {

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

# OLD VERSION
# stats_qPCR <- function(path_to_qPCR_res,
#                        out_dir_path,
#                        controls = c("18S", "GAPDH", "RPL13A"),
#                        stats_type = "ttest",
#                        stats_value = "log2_norm_exp",
#                        paired = F,
#                        var.equal = TRUE,
#                        alternative = "two.sided",
#                        stats_group = NULL,
#                        file_prefix = "") {
#
#   if (!stats_type %in% c("ttest", "anova")) {
#     stop("stats_type must be either 'ttest' or 'anova'.")
#   }
#
#   if (!stats_value %in% c("log2_norm_exp", "norm_exp")) {
#     stop("stats_value must be either 'log2_norm_exp' or 'norm_exp'.")
#   }
#
#   # Load dependent libraries----------------------------------------------------
#   suppressPackageStartupMessages(require(readr))
#   suppressPackageStartupMessages(require(dplyr))
#   suppressPackageStartupMessages(require(purrr))
#   suppressPackageStartupMessages(require(broom))
#   suppressPackageStartupMessages(require(rstatix))
#
#   # make outdir folder
#
#   dir.create(file.path(out_dir_path),
#              recursive = TRUE,
#              showWarnings = FALSE)
#
#   # Load a qPCR data containing log2 normalized expression values.
#   qpcr_data <- readr::read_tsv(path_to_qPCR_res,
#                                col_types = cols())
#
#   # Define factor for ordering the group label.
#   group_fct <- qpcr_data %>%
#     distinct(Group, Group_order) %>%
#     arrange(Group_order) %>%
#     pull(Group)
#
#   qpcr_data <- qpcr_data %>%
#     dplyr::mutate(Group = factor(Group, level = group_fct))
#
#   genes <- qpcr_data %>%
#     dplyr::filter(!Target %in% controls) %>%
#     dplyr::pull(Target) %>%
#     unique()
#
#   merged_stat_test <- NULL
#
#   if (stats_type == "ttest") {
#
#     for (i in genes) {
#
#       cat(paste0("Perform ttest on ",i,".\n"))
#
#       df <- qpcr_data %>%
#         dplyr::filter(Target == i)
#
#       if (1 %in% df$n | 0 %in% df$n) {
#
#         cat(paste0("At least one group has just 1 or 0 sample in ",i,".\n"))
#         cat(paste0("Cannot perform ttest on ",i,".\n"))
#
#       } else {
#
#         if (stats_value == "log2_norm_exp"){
#           stat_test <- t_test(log2_norm_exp ~ Group,
#                               paired = paired,
#                               var.equal = var.equal,
#                               alternative = alternative,
#                               data = df)
#
#           dummy_test <- t_test(norm_exp ~ Group, data = df) %>%
#             add_y_position() %>%
#             dplyr::select(group1,group2,y.position)
#
#           stat_test <- stat_test %>%
#             left_join(dummy_test, by = c("group1", "group2"))
#
#         } else if (stats_value == "norm_exp") {
#
#           stat_test <- t_test(norm_exp ~ Group,
#                               paired = paired,
#                               var.equal = var.equal,
#                               alternative = alternative,
#                               data = df) %>%
#             add_y_position()
#
#         }
#
#         stat_test <- stat_test %>% dplyr::mutate(Target = i, .before = group1)
#
#         merged_stat_test <- bind_rows(merged_stat_test, stat_test) %>%
#           dplyr::mutate(p.signif = case_when(p > 0.05 ~ "ns",
#                                              p <= 0.05 & p > 0.01 ~ "*",
#                                              p <= 0.01 & p > 0.001 ~ "**",
#                                              p <= 0.001 ~ "***"),
#                         .before = y.position)
#       }
#
#     }
#
#     cat("Calculation has been done.\n")
#     readr::write_tsv(merged_stat_test,
#                      file = paste0(out_dir_path, "/",
#                                    file_prefix,
#                                    "qPCR_ttest_result.tsv"))
#     cat("The table has been saved.\n")
#
#   } else if (stats_type == "anova") {
#
#     for (i in genes) {
#
#       cat(paste0("Perform anova on ",i,".\n"))
#
#       df <- qpcr_data %>%
#         dplyr::filter(Target == i)
#
#       if (1 %in% df$n | 0 %in% df$n) {
#
#         cat(paste0("At least one group has just 1 or 0 sample in ",i,".\n"))
#         cat(paste0("Cannot perform anova on ",i,".\n"))
#
#       } else {
#
#         if (stats_value == "log2_norm_exp") {
#           stat_test <- aov(log2_norm_exp ~ Group, data = df) %>%
#             tukey_hsd()
#
#           dummy_test <- t_test(norm_exp ~ Group, data = df) %>%
#             add_y_position() %>%
#             dplyr::select(group1,group2,y.position)
#
#           stat_test <- stat_test %>%
#             left_join(dummy_test, by = c("group1", "group2"))
#         } else if (stats_value == "norm_exp") {
#           stat_test <- aov(norm_exp ~ Group, data = df) %>%
#             tukey_hsd()
#
#           dummy_test <- t_test(norm_exp ~ Group, data = df) %>%
#             add_y_position() %>%
#             dplyr::select(group1,group2,y.position)
#
#           stat_test <- stat_test %>%
#             left_join(dummy_test, by = c("group1", "group2"))
#
#         }
#
#         stat_test <- stat_test %>% dplyr::mutate(Target = i, .before = group1)
#
#         merged_stat_test <- bind_rows(merged_stat_test, stat_test)
#
#       }
#     }
#     cat("Calculation has been done.\n")
#     readr::write_tsv(merged_stat_test,
#                      file = paste0(out_dir_path, "/",
#                                    file_prefix,
#                                    "qPCR_anova_tukey_result.tsv"))
#     cat("The table has been saved.\n")
#
#   }
#
# }
