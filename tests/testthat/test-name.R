test_that("named palettes are matched to group names", {
  palette <- c(
    Treatment_B = "#4DBBD5",
    Control = "#FFFFFF",
    Treatment_A = "#E64B35"
  )

  result <- molecularBio:::resolve_palette(
    palette,
    groups = c("Control", "Treatment_A", "Treatment_B")
  )

  expect_identical(
    result,
    c(
      Control = "#FFFFFF",
      Treatment_A = "#E64B35",
      Treatment_B = "#4DBBD5"
    )
  )
})

test_that("unnamed palettes remain positional", {
  result <- molecularBio:::resolve_palette(
    c("#FFFFFF", "#E64B35"),
    groups = c("Control", "Treatment")
  )

  expect_identical(
    result,
    c(Control = "#FFFFFF", Treatment = "#E64B35")
  )
})

test_that("invalid named palettes fail clearly", {
  expect_error(
    molecularBio:::resolve_palette(
      c(Control = "#FFFFFF"),
      groups = c("Control", "Treatment")
    ),
    "does not define a color for every group"
  )

  expect_error(
    molecularBio:::resolve_palette(
      c(Control = "not-a-color"),
      groups = "Control"
    ),
    "invalid R color"
  )
})

test_that("individual plot jitter is configurable and reproducible", {
  object <- structure(
    list(
      data = data.frame(
        Group = rep(c("Control", "Treatment"), each = 2),
        Target = "GeneA",
        norm_exp = c(0.9, 1.1, 1.8, 2.2),
        group_exp = rep(c(1, 2), each = 2),
        upper_error = rep(c(1.1, 2.2), each = 2),
        lower_error = rep(c(0.9, 1.8), each = 2),
        n = 2
      ),
      controls = "GAPDH",
      group_levels = c("Control", "Treatment"),
      reference_group = "Control"
    ),
    class = "qpcr_analysis"
  )
  colors <- c(Control = "#FFFFFF", Treatment = "#E64B35")

  jittered <- molecularBio:::build_gene_plot(
    object,
    gene = "GeneA",
    colors = colors,
    jitter = TRUE,
    jitter_seed = 1L
  )
  centered <- molecularBio:::build_gene_plot(
    object,
    gene = "GeneA",
    colors = colors,
    jitter = FALSE
  )

  expect_s3_class(jittered$layers[[3]]$position, "PositionJitter")
  expect_identical(jittered$layers[[3]]$position$seed, 1L)
  expect_s3_class(centered$layers[[3]]$position, "PositionIdentity")
})

test_that("technical replicates are detected and averaged", {
  raw_data <- data.frame(
    Sample = c("Sample1", "Sample1", "Sample1"),
    Target = "GeneA",
    Cq = c("20.0", "20.2", "Undetermined"),
    `Amp Status` = c("Amp", "Amp", "No Amp"),
    check.names = FALSE
  )

  expect_message(
    suppressWarnings(
      molecularBio:::collapse_technical_replicates(
        raw_data,
        technical_sd_threshold = 0.1
      )
    ),
    "Averaging technical replicates"
  )
  expect_warning(
    result <- suppressMessages(
      molecularBio:::collapse_technical_replicates(
        raw_data,
        technical_sd_threshold = 0.1
      )
    ),
    "SD exceeded"
  )

  expect_equal(result$Cq, 20.1)
  expect_equal(result$technical_n, 3L)
  expect_equal(result$technical_valid_n, 2L)
  expect_equal(result$technical_not_determined_n, 1L)
  expect_equal(result$technical_sd, stats::sd(c(20, 20.2)))
})

test_that("absence of technical replicates is reported", {
  raw_data <- data.frame(
    Sample = c("Sample1", "Sample2"),
    Target = "GeneA",
    Cq = c(20, 21),
    `Amp Status` = "Amp",
    check.names = FALSE
  )

  expect_message(
    result <- molecularBio:::collapse_technical_replicates(raw_data),
    "No technical replicates were provided"
  )
  expect_equal(result$technical_n, c(1L, 1L))
})

test_that("technical replicate error mode rejects repeated pairs", {
  raw_data <- data.frame(
    Sample = c("Sample1", "Sample1"),
    Target = "GeneA",
    Cq = c(20, 20.1),
    `Amp Status` = "Amp",
    check.names = FALSE
  )

  expect_error(
    molecularBio:::collapse_technical_replicates(
      raw_data,
      technical_replicates = "error"
    ),
    "Technical replicates were detected"
  )
})
