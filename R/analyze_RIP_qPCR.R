
analyze_binding_qPCR <- function(
    path,
    out_dir = ".",
    plot = c("ip_input", "fold"),
    max_ct = 40,
    dpi = 300,
    width = 8,
    height = 4
) {

  plot <- match.arg(plot)

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
  })

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dat <- read_csv(
    path,
    show_col_types = FALSE
  ) %>%
    mutate(
      Cq = trimws(as.character(Cq)),
      Cq = case_when(
        is.na(Cq) ~ as.character(max_ct),
        Cq == "" ~ as.character(max_ct),
        grepl("^Undetermined$", Cq, ignore.case = TRUE) ~ as.character(max_ct),
        TRUE ~ Cq
      ),
      Cq = readr::parse_number(Cq)
    )

  res <- dat %>%
    pivot_wider(
      id_cols = c(Target_order, Target, `Input %`),
      names_from = Type,
      values_from = Cq
    ) %>%
    mutate(
      dilution_factor = 100 / `Input %`,
      adjusted_input = Input + log2(dilution_factor),
      IP_Input = 2^(-(IP - adjusted_input)),
      IgG_Input = 2^(-(IgG - adjusted_input)),
      Fold_enrichment = IP_Input / IgG_Input,
      IgG_ND = IgG >= max_ct
    ) %>%
    arrange(Target_order)

  target_levels <- res %>%
    distinct(Target_order, Target) %>%
    arrange(Target_order) %>%
    pull(Target)

  res <- res %>%
    mutate(
      Target = factor(Target, levels = target_levels)
    )

  write_tsv(
    res,
    file.path(
      out_dir,
      "binding_qPCR_results.tsv"
    )
  )

  if (plot == "fold") {

    g <- ggplot(
      res,
      aes(
        x = Target,
        y = Fold_enrichment
      )
    ) +
      geom_col(
        fill = "#4C72B0",
        colour = "black",
        width = 0.7
      ) +
      geom_text(
        data = filter(res, IgG_ND),
        aes(label = "ND"),
        colour = "red",
        vjust = -0.4,
        size = 4
      ) +
      theme_classic(base_size = 14) +
      labs(
        x = NULL,
        y = "Fold enrichment"
      ) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )
      )

  } else {

    plot_df <- res %>%
      select(
        Target,
        IP_Input,
        IgG_Input
      ) %>%
      pivot_longer(
        cols = c(IP_Input, IgG_Input),
        names_to = "Fraction",
        values_to = "Value"
      )

    g <- ggplot(
      plot_df,
      aes(
        x = Target,
        y = Value,
        fill = Fraction
      )
    ) +
      geom_col(
        position = position_dodge(0.8),
        width = 0.7,
        colour = "black"
      ) +
      scale_fill_manual(
        values = c(
          IP_Input = "#4C72B0",
          IgG_Input = "grey75"
        )
      ) +
      theme_classic(base_size = 14) +
      labs(
        x = NULL,
        y = "Fraction of input"
      ) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )
      )

  }

  ggsave(
    filename = file.path(
      out_dir,
      paste0(plot, ".svg")
    ),
    plot = g,
    width = width,
    height = height
  )

  ggsave(
    filename = file.path(
      out_dir,
      paste0(plot, ".tiff")
    ),
    plot = g,
    width = width,
    height = height,
    dpi = dpi
  )

  invisible(res)

}

plot_binding_normalized <- function(
    path,
    out_dir = ".",
    max_ct = 40,
    palette = c(
      Input = "grey90",
      IP = "#3B7DDD",
      IgG = "grey50"
    ),
    width = 3,
    height = 4,
    dpi = 300
) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
    library(purrr)
  })

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dat <- read_csv(
    path,
    show_col_types = FALSE
  ) %>%
    mutate(
      Cq = trimws(as.character(Cq)),
      Cq = case_when(
        is.na(Cq) ~ as.character(max_ct),
        Cq == "" ~ as.character(max_ct),
        grepl("^Undetermined$", Cq, ignore.case = TRUE) ~ as.character(max_ct),
        TRUE ~ Cq
      ),
      Cq = readr::parse_number(Cq)
    )

  res <- dat %>%
    pivot_wider(
      id_cols = c(Target, `Input %`),
      names_from = Type,
      values_from = Cq
    ) %>%
    mutate(
      dilution_factor = 100 / `Input %`,
      adjusted_input = Input + log2(dilution_factor),
      Input_Input = 1,
      IP_Input = 2^(-(IP - adjusted_input)),
      IgG_Input = 2^(-(IgG - adjusted_input)),
      IgG_ND = IgG >= max_ct
    )

  walk(
    unique(res$Target),
    function(gene) {

      df <- res %>%
        filter(Target == gene) %>%
        select(
          Input_Input,
          IP_Input,
          IgG_Input,
          IgG_ND
        ) %>%
        pivot_longer(
          cols = c(
            Input_Input,
            IP_Input,
            IgG_Input
          ),
          names_to = "Fraction",
          values_to = "Value"
        ) %>%
        mutate(
          Fraction = factor(
            Fraction,
            levels = c(
              "Input_Input",
              "IP_Input",
              "IgG_Input"
            ),
            labels = c(
              "Input",
              "IP",
              "IgG"
            )
          )
        )

      nd_flag <- res %>%
        filter(Target == gene) %>%
        pull(IgG_ND)

      df<- df %>%
        mutate(
          Value = if_else(
            Fraction == "IgG" & IgG_ND,
            NA_real_,
            Value
          )
        )

      p <- ggplot(
        df,
        aes(
          Fraction,
          Value,
          fill = Fraction
        )
      ) +
        geom_col(
          width = 0.7,
          colour = "black"
        ) +
        scale_fill_manual(
          values = palette
        ) +
        theme_classic(base_size = 14) +
        labs(
          title = gene,
          x = NULL,
          y = "Relative enrichment\n(Input = 1)"
        ) +
        theme(
          legend.position = "none",
          plot.title = element_text(
            hjust = 0.5
          )
        )

      if (nd_flag) {

        p <- p +
          annotate(
            "text",
            x = 3,
            y = max(df$Value, na.rm = TRUE) * 1.05,
            label = "ND",
            colour = "red",
            fontface = "bold",
            size = 5
          )

      }

      ggsave(
        filename = file.path(
          out_dir,
          paste0(gene, "_normalized.svg")
        ),
        plot = p,
        width = width,
        height = height
      )

      ggsave(
        filename = file.path(
          out_dir,
          paste0(gene, "_normalized.tiff")
        ),
        plot = p,
        width = width,
        height = height,
        dpi = dpi
      )

    }
  )

  invisible(res)

}


plot_binding_percent_input <- function(
    path,
    out_dir = ".",
    max_ct = 40,
    palette = c(
      Input = "grey90",
      IP = "#3B7DDD",
      IgG = "grey50"
    ),
    width = 3,
    height = 4,
    dpi = 300
) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
    library(purrr)
  })

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dat <- read_csv(
    path,
    show_col_types = FALSE
  ) %>%
    mutate(
      Cq = trimws(as.character(Cq)),
      Cq = case_when(
        is.na(Cq) ~ as.character(max_ct),
        Cq == "" ~ as.character(max_ct),
        grepl("^Undetermined$", Cq, ignore.case = TRUE) ~
          as.character(max_ct),
        TRUE ~ Cq
      ),
      Cq = parse_number(Cq)
    )

  res <- dat %>%
    group_by(Target) %>%
    summarise(
      Input = Cq[Type == "Input"][1],
      IP = Cq[Type == "IP"][1],
      IgG = Cq[Type == "IgG"][1],
      `Input %` = first(`Input %`),
      .groups = "drop"
    ) %>%
    mutate(

      dilution_factor = 100 / `Input %`,

      adjusted_input = Input + log2(dilution_factor),

      Input_percent = dilution_factor,

      IP_percent =
        2^(-(IP - adjusted_input)),

      IgG_percent =
        2^(-(IgG - adjusted_input)),

      IgG_ND = IgG >= max_ct
    )

  walk(
    seq_len(nrow(res)),
    function(i) {

      gene <- res$Target[i]

      df <- tibble(
        Fraction = factor(
          c("Input", "IP", "IgG"),
          levels = c("Input", "IP", "IgG")
        ),
        Value = c(
          res$Input_percent[i],
          res$IP_percent[i],
          res$IgG_percent[i]
        )
      )

      if (res$IgG_ND[i]) {
        df$Value[df$Fraction == "IgG"] <- NA
      }

      ymax <- max(df$Value, na.rm = TRUE)

      p <-
        ggplot(
          df,
          aes(
            Fraction,
            Value,
            fill = Fraction
          )
        ) +
        geom_col(
          colour = "black",
          width = 0.7,
          na.rm = TRUE
        ) +
        scale_fill_manual(
          values = palette
        ) +
        labs(
          title = gene,
          x = NULL,
          y = "% Input"
        ) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "none",
          plot.title = element_text(
            hjust = 0.5
          )
        )

      if (res$IgG_ND[i]) {

        p <-
          p +
          annotate(
            "text",
            x = 3,
            y = ymax * 1.05,
            label = "ND",
            colour = "red",
            fontface = "bold",
            size = 5
          )

      }

      ggsave(
        file.path(
          out_dir,
          paste0(gene, "_percent_input.svg")
        ),
        p,
        width = width,
        height = height
      )

      ggsave(
        file.path(
          out_dir,
          paste0(gene, "_percent_input.tiff")
        ),
        p,
        width = width,
        height = height,
        dpi = dpi
      )

    }

  )

  invisible(res)

}
