#' Analyze LucAssay data
#'
#' @description
#' `analyze_LucAssay()`
#'
#'
#'@param path_to_LucAssay_res path to a csv file produced by QuantStudio.
#'@param path_to_group_label path to a csv file including meta data.
#'@param out_dir_path path to an output directory.
#'@param save_final_table A boolean.
#'@param file_prefix prefix for files.
#'@export

analyze_LucAssay <- function(path_to_firefly_res,
                             path_to_renilla_res,
                         path_to_group_label,
                         out_dir_path,
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

  # Load a raw LucAssay data containing signals
  raw_data_f <- readr::read_csv(path_to_firefly_res,
                              comment = "",
                              skip = 0,
                              col_types = cols())
  raw_data_r <- readr::read_csv(path_to_renilla_res,
                              comment = "",
                              skip = 0,
                              col_types = cols())

  # Load a file containing group label and group order info.
  group_label <- readr::read_csv(path_to_group_label,
                                 comment = "",
                                 skip = 0,
                                 col_types = cols())


  # Calculate mean Firefly_Luc value.
  summarize_data <- raw_data_f %>%
    dplyr::left_join(raw_data_r, by = "Sample") %>%
    dplyr::left_join(group_label[1:2], by = "Sample") %>%
    dplyr::relocate(Group, .after = "Sample") %>%
    dplyr::group_by(Group) %>%
    dplyr::mutate(
      FR_Ratio = Firefly_Luc_Signal/Renilla_Luc_Signal,
      Ave_Firefly_Luc_Signal = mean(Firefly_Luc_Signal, na.rm = TRUE),
      Ave_Renilla_Luc_Signal = mean(Renilla_Luc_Signal, na.rm = TRUE),
      Ave_FR_Ratio = mean(FR_Ratio, na.rm = TRUE),
      n = sum(!is.na(Sample)),
      SD_FL = if_else(n == 1, 0, sd(Firefly_Luc_Signal, na.rm = TRUE)),
      SEM_FL = SD_FL/sqrt(n),
      SD_RL = if_else(n == 1, 0, sd(Renilla_Luc_Signal, na.rm = TRUE)),
      SEM_RL = SD_RL/sqrt(n),
      SD_FR_Ratio = if_else(n == 1, 0, sd(FR_Ratio, na.rm = TRUE)),
      SEM_FR_Ratio = SD_FR_Ratio/sqrt(n))

  # Define factor for ordering the group label.
  group_fct <- group_label %>%
    distinct(Group, Group_order) %>%
    arrange(Group_order) %>%
    pull(Group)

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

  # color scheme
  # 1 = white - grey
  # 2 = red
  # 3 = blue

  # Define functions to create color gradients.
  colfunc_grey <- colorRampPalette(c("white", "grey80"))
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
    grey_palette <- "white"
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

  plot_data <- summarize_data %>%
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
                            "analyzed_LucAssay_res.tsv")
              )
  }
}
