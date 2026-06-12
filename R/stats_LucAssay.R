#' Perform statistical analysis on LucAssay data
#'
#' @description
#' `stats_LucAssay()` analyzes a tsv file produced by analyze_LucAssay(). \n
#'
#'@param path_to_LucAssay_res path to a tsv file produced by analyze_LucAssay().
#'@param out_dir_path path to an output directory.
#'@param file_prefix prefix for files.
#'@export

stats_LucAssay <- function(path_to_LucAssay_res,
                       out_dir_path,
                       stats_type = "ttest",
                       paired = F,
                       var.equal = TRUE,
                       alternative = "two.sided",
                       stats_group = NULL,
                       file_prefix = "") {

  if (!stats_type %in% c("ttest", "anova")) {
    stop("stats_type must be either 'ttest' or 'anova'.")
  }

  # Load dependent libraries----------------------------------------------------
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(tidyr))
  suppressPackageStartupMessages(require(purrr))
  suppressPackageStartupMessages(require(broom))
  suppressPackageStartupMessages(require(rstatix))

  # make outdir folder

  dir.create(file.path(out_dir_path),
             recursive = TRUE,
             showWarnings = FALSE)

  # Load a LucAssay data.
  LucAssay_data <- readr::read_tsv(path_to_LucAssay_res,
                               col_types = cols())

  # Define factor for ordering the group label.
  group_fct <- LucAssay_data %>%
    distinct(Group, Group_order) %>%
    arrange(Group_order) %>%
    pull(Group)

  LucAssay_data <- LucAssay_data %>%
    dplyr::mutate(Group = factor(Group, level = group_fct))

  merged_stat_test <- NULL

  records <- c("FireFly_Luc_Signal", "Renilla_Luc_Signal", "FR_Ratio")


  if (stats_type == "ttest") {

    cat(paste0("Run ttest.\n"))

    if (1 %in% LucAssay_data$n | 0 %in% LucAssay_data$n) {

      cat(paste0("At least one group has just 1 or 0 sample.\n"))
      cat(paste0("Cannot perform ttest.\n"))

    } else {

      stat_test <- LucAssay_data %>%
        dplyr::select(Group,
                      Firefly_Luc_Signal,
                      Renilla_Luc_Signal,
                      FR_Ratio) %>%
        tidyr::pivot_longer(-Group,
                            names_to = "variables",
                            values_to = "value") %>%
        tidyr::nest(data = -variables) %>%
        dplyr::mutate(
          ttest_results = purrr::map(data, ~t_test(value ~ Group,
                                                   data = .x,
                                                   paired = paired,
                                                   var.equal = var.equal,
                                                   alternative = alternative)),
          ttest_results = purrr::map2(ttest_results, data,
                                      ~ .x %>%
                                        add_y_position(data = .y,
                                                       formula = value ~ Group))
          ) %>%
        select(variables, ttest_results) %>%
        unnest(cols = c(ttest_results))

    }

    cat("Calculation has been done.\n")
    readr::write_tsv(stat_test,
                     file = paste0(out_dir_path, "/",
                                   file_prefix,
                                   "LucAssay_ttest_result.tsv"))
    cat("The table has been saved.\n")

  } else if (stats_type == "anova") {

    if (1 %in% LucAssay_data$n | 0 %in% LucAssay_data$n) {

      cat(paste0("At least one group has just 1 or 0 sample.\n"))
      cat(paste0("Cannot perform anova.\n"))

    } else {
      cat(paste0("Run anova.\n"))
      stat_test <- LucAssay_data %>%
        dplyr::select(Group,
                      Firefly_Luc_Signal,
                      Renilla_Luc_Signal,
                      FR_Ratio) %>%
        tidyr::pivot_longer(-Group,
                            names_to = "variables",
                            values_to = "value") %>%
        tidyr::nest(data = -variables) %>%
        dplyr::mutate(anova_res = purrr::map(data, ~aov(value ~ Group, .)),
                      tukey_res = purrr::map(anova_res, ~tukey_hsd(.))
                      ) %>%
        select(variables, tukey_res) %>%
        unnest(cols = c(tukey_res))

      dummy_test <- LucAssay_data %>%
        dplyr::select(Group,
                      Firefly_Luc_Signal,
                      Renilla_Luc_Signal,
                      FR_Ratio) %>%
        tidyr::pivot_longer(-Group,
                            names_to = "variables",
                            values_to = "value") %>%
        tidyr::nest(data = -variables) %>%
        dplyr::mutate(
          ttest_results = purrr::map(data, ~t_test(value ~ Group,
                                                   data = .x,
                                                   paired = paired,
                                                   var.equal = var.equal,
                                                   alternative = alternative)),
          ttest_results = purrr::map2(ttest_results, data,
                                      ~ .x %>%
                                        add_y_position(data = .y,
                                                       formula = value ~ Group))
        ) %>%
        select(variables, ttest_results) %>%
        unnest(cols = c(ttest_results)) %>%
        dplyr::select(variables, group1, group2, y.position)

      stat_test <- stat_test %>%
        left_join(dummy_test, by = c("variables", "group1", "group2"))

    }

    cat("Calculation has been done.\n")
    readr::write_tsv(stat_test,
                     file = paste0(out_dir_path, "/",
                                   file_prefix,
                                   "LucAssay_anova_tukey_result.tsv"))
    cat("The table has been saved.\n")

  }

}
