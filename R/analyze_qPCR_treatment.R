#' Analyze qPCR data
#'
#' @description
#' `analyze_qPCR_treatment()` analyzes a csv file produced by QuantStudio based on
#' the ddCT method. This function implements the method described in
#' this reference.
#' Reference https://www.cell.com/trends/biotechnology/fulltext/S0167-7799(18)30342-1
#'
#'
#'@param path_to_qPCR_res path to a csv file produced by QuantStudio.
#'@param path_to_group_label path to a csv file including meta data.
#'@param out_dir_path path to an output directory.
#'@param controls list of internal controls.
#'      c("18S", "GAPDH", "RPL13A") (default)
#'@param save_final_table A boolean.
#'@param file_prefix prefix for files.
#'@export

analyze_qPCR_treatment <- function(path_to_qPCR_res,
                         path_to_group_label,
                         out_dir_path,
                         controls = c("18S", "GAPDH", "RPL13A"),
                         save_final_table = TRUE,
                         file_prefix = "") {

  # Load dependent libraries----------------------------------------------------
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(ggtext))

  # Load and modify data for plotting ------------------------------------------

  # make outdir folder

  dir.create(file.path(out_dir_path),
             recursive = TRUE,
             showWarnings = FALSE)

  # Load a raw qPCR data containing Ct values.
  raw_data <- readr::read_csv(path_to_qPCR_res,
                              comment = "",
                              skip = 0,
                              col_types = cols())

  # Load a file containing group label and group order info.
  group_label <- readr::read_csv(path_to_group_label,
                                 comment = "",
                                 skip = 0,
                                 col_types = cols())

  # Calculate mean Ct value.
  summarize_data <- raw_data %>%
    dplyr::left_join(group_label[1:3], by = "Sample") %>%
    dplyr::select(Sample, Group, treatment_conc, Target, Cq, `Amp Status`) %>%
    dplyr::group_by(Group, Target, treatment_conc) %>%
    dplyr::mutate(Cq = dplyr::case_when(`Amp Status` == "Amp" ~ as.numeric(Cq)),
                  Cq_mean = mean(Cq, na.rm = TRUE))

  # Extract mean Ct values of the control group.
  base_mean <- summarize_data %>%
    dplyr::filter(treatment_conc == 0) %>%
    dplyr::distinct(Group, Target, Cq_mean)

  # By using base_mean, calculate dCt('Average Ct of control' - 'sample Ct').
  summarize_data_dCT <- summarize_data %>%
    dplyr::left_join(base_mean, by = c("Group", "Target"), suffix = c("", "_base")) %>%
    dplyr::mutate(dCt = Cq_mean_base - Cq,
                  RQ = 2**dCt)

  # By using RQ of internal control genes, calculate norm factor.
  norm_factors <- summarize_data_dCT %>%
    dplyr::filter(Target %in% controls) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(norm_factor = exp(mean(log(RQ), na.rm = TRUE))) %>%
    dplyr::select(Sample, norm_factor) %>%
    dplyr::distinct(Sample, norm_factor)

  # By using norm factor, calculate normalized RQ, mean of RQ, SD, SEM, etc.
  final_data <- summarize_data_dCT %>%
    dplyr::left_join(norm_factors, by = "Sample") %>%
    dplyr::mutate(norm_exp = RQ/norm_factor) %>%
    dplyr::mutate(log2_norm_exp = log2(norm_exp)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Group, Target, treatment_conc) %>%
    dplyr::mutate(group_exp = exp(mean(log(norm_exp), na.rm = TRUE)),
                  log2_group_exp = log2(group_exp),
                  Not_determined_n = sum(is.na(Cq)),
                  n = sum(!is.na(Cq)),
                  total_n = Not_determined_n + n,
                  SD = if_else(n == 1, 0, sd(log2_norm_exp, na.rm = TRUE)),
                  SEM = SD/sqrt(n),
                  upper_error = 2**(log2_group_exp + SEM),
                  lower_error = 2**(log2_group_exp - SEM))

  # Define factor for ordering the group label.
  group_fct <- group_label %>%
    distinct(Group, Group_order) %>%
    arrange(Group_order) %>%
    pull(Group)


  # color scheme
  # 1 = white - grey
  # 2 = red
  # 3 = blue

  # Define functions to create color gradients.
  colfunc_grey <- colorRampPalette(c("black", "grey80"))
  colfunc_red <- colorRampPalette(c("red", "pink"))
  colfunc_blue <- colorRampPalette(c("blue", "lightblue"))

  # Calculate number of groups within each color group, and rank them.
  color_scheme <- group_label %>%
    dplyr::distinct(Group,Group_order, color_scheme) %>%
    dplyr::arrange(Group_order) %>%
    dplyr::group_by(color_scheme) %>%
    dplyr::mutate(color_n=n()) %>%
    dplyr::mutate(color_order = row_number())

  # Make color palette

  if (1 %in% color_scheme$color_scheme){
    grey_palette <- colfunc_grey(
      color_scheme %>%
        filter(color_scheme == 1) %>%
        pull(color_n) %>%
        unique()
    )
  } else {
    grey_palette <- "black"
  }

  if (2 %in% color_scheme$color_scheme){
    red_palette <- colfunc_red(
      color_scheme %>%
        filter(color_scheme == 2) %>%
        pull(color_n) %>%
        unique()
    )
  } else {
    red_palette <- "red"
  }

  if (3 %in% color_scheme$color_scheme){
    blue_palette <- colfunc_blue(
      color_scheme %>%
        filter(color_scheme == 3) %>%
        pull(color_n) %>%
        unique()
    )
  } else {
    blue_palette <- "blue"
  }

  plot_data <- final_data %>%
    dplyr::left_join(color_scheme, by = "Group") %>%
    dplyr::mutate(Group = factor(Group, level = group_fct)) %>%
    dplyr::arrange(Group) %>%
    dplyr::mutate(color = case_when(color_scheme == 1 ~ grey_palette[color_order],
                                    color_scheme == 2 ~ red_palette[color_order],
                                    color_scheme == 3 ~ blue_palette[color_order],
                                    TRUE ~ "black"))

  # Save table------------------------------------------------------------------

  if (save_final_table == TRUE) {
    write_tsv(plot_data,
              file = paste0(out_dir_path, "/",
                            file_prefix,
                            "analyzed_qPCR_res.tsv")
              )
  }
}
