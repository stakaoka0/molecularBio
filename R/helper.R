# helpers.R

#' Geometric mean
#'
#' @keywords internal
geo_mean <- function(
    x,
    na.rm = TRUE
) {
  exp(
    mean(
      log(x),
      na.rm = na.rm
    )
  )
}

#' Run statistics for one gene
#'
#' @keywords internal
#'

run_gene_statistics <- function(
  df,
  formula,
  stats_type = c("ttest", "anova"),
  paired = FALSE,
  var.equal = TRUE,
  alternative = "two.sided"
) {
  stats_type <- match.arg(
    stats_type
  )
  gene <- unique(df$Target)
  group_sizes <- df |>
    dplyr::distinct(
      Group,
      n
    )
  if (any(group_sizes$n < 2)) {
    cli::cli_warn(
      paste0(
        "Skipping ",
        gene,
        ": at least one group has < 2 samples."
      )
    )
    return(
      tibble::tibble()
    )
  }

  if (stats_type == "ttest") {

    res <- rstatix::t_test(
      formula,
      data = df,
      paired = paired,
      var.equal = var.equal,
      alternative = alternative
    )
    ypos <- compute_y_positions(
      df
    )
    res <- dplyr::left_join(
      res,
      ypos,
      by = c(
        "group1",
        "group2"
      )
    )
    res <- res |>
      dplyr::mutate(
        Target = gene,
        p.signif =
          get_significance_label(p),
        .before = 1
      )
  } else {
    res <- stats::aov(
      formula,
      data = df
    ) |>
      rstatix::tukey_hsd()
    ypos <- compute_y_positions(
      df
    )
    res <- dplyr::left_join(
      res,
      ypos,
      by = c(
        "group1",
        "group2"
      )
    )
    res <- res |>
      dplyr::mutate(
        Target = gene,
        .before = 1
      )
  }
  res
}

#' Create significance labels
#'
#' @keywords internal

get_significance_label <- function(p) {
  dplyr::case_when(
    p > 0.05 ~ "ns",
    p > 0.01 ~ "*",
    p > 0.001 ~ "**",
    TRUE ~ "***"
  )
}


filter_stats_results <- function(
    stat_df,
    comparison = "all",
    reference_group = NULL,
    comparisons = NULL
) {
  if (comparison == "all") {
    return(stat_df)
  }
  if (comparison == "significant") {
    p_col <- if ("p.adj" %in% names(stat_df)) {
      "p.adj"
    } else {
      "p"
    }
    return(
      stat_df %>%
        dplyr::filter(
          .data[[p_col]] <= 0.05
        )
    )
  }
  if (comparison == "control") {
    return(
      stat_df %>%
        dplyr::filter(
          group1 == reference_group |
            group2 == reference_group
        )
    )
  }
  if (comparison == "custom") {
    keep <- purrr::map_dfr(
      comparisons,
      ~ tibble::tibble(
        group1 = .x[1],
        group2 = .x[2]
      )
    )
    return(
      dplyr::semi_join(
        stat_df,
        keep,
        by = c("group1", "group2")
      )
    )
  }

}


#' Generate y-axis label
#'
#' @keywords internal
make_y_label <- function(
    controls,
    show_controls = TRUE
) {
  if (!show_controls) {
    return("Relative Expression")
  }
  paste0(
    "Relative Expression\n",
    "(Normalized by ",
    paste(controls, collapse = ", "),
    ")"
  )
}



#' Build grouped color palette
#'
#' @keywords internal
#'

build_grouped_palette <- function(
    groups,
    group_info
) {

  family_colors <- list(
    "1" = c("white", "grey70"),
    "2" = c("#D73027", "#F4A582"),
    "3" = c("#4575B4", "#92C5DE"),
    "4" = c("#1A9850", "#A6D96A"),
    "5" = c("#7B3294", "#C2A5CF"),
    "6" = c("#E66101", "#FDB863"),
    "7" = c("#8C510A", "#D8B365"),
    "8" = c("#008080", "#80CDC1")
  )
  palette_tbl <- group_info |>
    dplyr::distinct(
      Group,
      Group_order,
      color_scheme
    ) |>
    dplyr::arrange(
      Group_order
    )
  palette_tbl <- palette_tbl |>
    dplyr::group_by(
      color_scheme
    ) |>
    dplyr::mutate(
      n_color = dplyr::n(),
      rank_color =
        dplyr::row_number()
    ) |>
    dplyr::ungroup()
  palette_tbl <- palette_tbl |>
    dplyr::rowwise() |>
    dplyr::mutate(
      color =
        grDevices::colorRampPalette(
          family_colors[[as.character(color_scheme)]]
        )(n_color)[rank_color]
      ) |>
    dplyr::ungroup()
  cols <- palette_tbl$color
  names(cols) <- palette_tbl$Group
  cols
}

#' Build khroma  palette
#'
#' @keywords internal
#'

build_khroma_palette <- function(
  groups,
  type = "bright"
) {
  cols <- khroma::colour(
    type
  )(
    length(groups)
  )
  names(cols) <- groups
  cols
}

#' Resolve color palette
#'
#' @keywords internal
#'

resolve_palette <- function(
    palette,
    groups,
    group_info = NULL
) {
  n_groups <- length(groups)
  if (
    is.character(palette) &&
    length(palette) == 1 &&
    (is.null(names(palette)) || !nzchar(names(palette)))
  ) {

    colors <- switch(
      palette,
      grouped =
        build_grouped_palette(
          groups = groups,
          group_info = group_info
        ),
      khroma_bright =
        build_khroma_palette(
          groups,
          "bright"
        ),
      khroma_muted =
        build_khroma_palette(
          groups,
          "muted"
        ),
      khroma_vibrant =
        build_khroma_palette(
          groups,
          "vibrant"
        ),
      khroma_light =
        build_khroma_palette(
          groups,
          "light"
        ),
      npg =
        ggsci::pal_npg("nrc")(
          n_groups
        ),
      aaas =
        ggsci::pal_aaas("default")(
          n_groups
        ),
      lancet =
        ggsci::pal_lancet()(n_groups),
      jama =
        ggsci::pal_jama()(n_groups),
      nejm =
        ggsci::pal_nejm()(n_groups),
      jco =
        ggsci::pal_jco()(n_groups),
      okabe =
        ggokabeito::palette_okabe_ito()[
          seq_len(n_groups)
        ],
      viridis =
        viridisLite::viridis(
          n_groups
        ),
      scales::hue_pal()(n_groups)
    )

    names(colors) <- groups
    return(colors)
  }
  palette_names <- names(palette)
  has_names <- !is.null(palette_names) &&
    any(nzchar(palette_names))
  if (has_names && any(!nzchar(palette_names))) {
    cli::cli_abort(
      "{.arg palette} must be either fully named or completely unnamed."
    )
  }
  if (has_names && anyDuplicated(palette_names)) {
    duplicated_names <- unique(
      palette_names[duplicated(palette_names)]
    )
    cli::cli_abort(
      c(
        "{.arg palette} contains duplicate group names.",
        "x Duplicated: {.val {duplicated_names}}"
      )
    )
  }
  if (has_names) {
    missing_groups <- setdiff(groups, palette_names)
    if (length(missing_groups) > 0) {
      cli::cli_abort(
        c(
          "{.arg palette} does not define a color for every group.",
          "x Missing: {.val {missing_groups}}"
        )
      )
    }
    colors <- palette[groups]
  } else {
    if (length(palette) < n_groups) {
      cli::cli_abort(
        "{.arg palette} contains fewer colors than groups."
      )
    }
    colors <- palette[seq_len(n_groups)]
    names(colors) <- groups
  }
  tryCatch(
    grDevices::col2rgb(unname(colors)),
    error = function(error) {
      cli::cli_abort(
        "{.arg palette} contains an invalid R color.",
        parent = error
      )
    }
  )
  colors
}


#' Compute y positions for stat annotations
#'
#' @keywords internal

compute_y_positions <- function(df) {
  rstatix::t_test(
    norm_exp ~ Group,
    data = df
  ) |>
    rstatix::add_y_position() |>
    dplyr::select(
      group1,
      group2,
      y.position
    )
}


#' Save plot helper
#'
#' @keywords internal

save_plot <- function(
  plot,
  filename_base,
  width,
  height,
  dpi,
  save_tiff = TRUE,
  save_svg = TRUE
) {

  if (save_tiff) {
    ggplot2::ggsave(
      filename = paste0(filename_base, ".tiff"),
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      compression = "lzw"
    )
  }
  if (save_svg) {
    ggplot2::ggsave(
      filename = paste0(filename_base, ".svg"),
      plot = plot,
      width = width,
      height = height,
      device = svglite::svglite
    )
  }

}
