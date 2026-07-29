
<!-- README.md is generated from README.Rmd. Please edit that file -->

# molecularBio

<!-- badges: start -->
<!-- badges: end -->

`molecularBio` provides a reproducible workflow for routine qPCR
analysis. It reads QuantStudio results, normalizes relative expression
against one or more reference genes, performs optional gene-wise
statistical tests, and creates publication-ready plots.

The calculation workflow is based on the relative-quantification
procedure described by Taylor *et al.* (2019), particularly the sequence
summarized in their Figure 5. See [Method correspondence and current
assumptions](#method-correspondence-and-current-assumptions) below.

The package currently focuses on standard qPCR experiments. Actinomycin
D, treatment-series, RIP-qPCR, and luciferase-assay workflows are not
part of the current `main` branch.

## Installation

Install the package from GitHub:

``` r
install.packages("remotes")
remotes::install_github("stakaoka0/molecularBio")
```

## Input files

### QuantStudio results

The qPCR CSV file must contain:

| Column       | Description                                |
|--------------|--------------------------------------------|
| `Sample`     | Sample identifier                          |
| `Target`     | Target or reference-gene name              |
| `Cq`         | Quantification cycle                       |
| `Amp Status` | Measurements equal to `"Amp"` are retained |

Values not marked `"Amp"` are treated as missing Cq measurements.
Repeated rows with the same `Sample` and `Target` are interpreted as
technical replicates and averaged before normalization. Consequently,
`Sample` must be a unique biological-sample identifier rather than a
group label.

### Sample metadata

The metadata CSV file must contain:

| Column         | Description                                                   |
|----------------|---------------------------------------------------------------|
| `Sample`       | Sample identifier matching the qPCR file                      |
| `Group`        | Experimental group                                            |
| `Group_order`  | Plotting order; the group numbered `1` is the reference group |
| `color_scheme` | Optional grouped-palette family                               |

Exactly one group must have `Group_order = 1`. If `color_scheme` is
omitted, the default grouped grey palette is used.

## Complete workflow

`qPCR()` runs analysis, optional statistics, and plotting in one call:

``` r
library(molecularBio)

result <- qPCR(
  path_to_qPCR_res = "data/quantstudio_results.csv",
  path_to_group_label = "data/sample_metadata.csv",
  out_dir_path = "results",
  controls = c("18S", "GAPDH", "RPL13A"),
  technical_replicates = "average",
  technical_sd_threshold = 0.2,
  stats_type = "anova",
  plot_all = TRUE,
  plot_individual = TRUE
)
```

Set `stats_type` to `"none"`, `"ttest"`, or `"anova"`. ANOVA is followed
by Tukey pairwise comparisons. Reference genes are excluded from
statistical testing.

The returned object contains the analyzed data and, when requested, the
statistical results:

``` r
result
summary(result)
result$data
result$stats
```

## Step-by-step workflow

The same workflow can be run in separate stages.

### Analyze expression

``` r
analysis <- analyze_qPCR(
  path_to_qPCR_res = "data/quantstudio_results.csv",
  path_to_group_label = "data/sample_metadata.csv",
  controls = c("18S", "GAPDH", "RPL13A")
)
```

For each target, the analysis:

1.  calculates the mean Cq in the reference group;
2.  calculates relative quantity as `2^(reference mean Cq - sample Cq)`;
3.  calculates a sample normalization factor from the geometric mean of
    the reference-gene relative quantities;
4.  reports normalized expression, group geometric means, SD, SEM, and
    error limits.

### Technical replicates

By default, repeated `Sample`-`Target` combinations are treated as
technical replicates:

``` r
analysis <- analyze_qPCR(
  path_to_qPCR_res = "data/quantstudio_results.csv",
  path_to_group_label = "data/sample_metadata.csv",
  technical_replicates = "average",
  technical_sd_threshold = 0.2
)
```

Only Cq values with `Amp Status == "Amp"` contribute to the technical
mean. The analyzed table retains:

- `technical_n`: total reactions for the sample-target pair;
- `technical_valid_n`: reactions contributing to the mean;
- `technical_not_determined_n`: excluded or undetermined reactions;
- `technical_sd`: SD among valid technical-replicate Cq values.

A warning identifies pairs whose technical SD exceeds
`technical_sd_threshold`. If no repeated pairs are present, the package
reports `"No technical replicates were provided."` To reject repeated
pairs instead of averaging them, use `technical_replicates = "error"`.

## Method correspondence and current assumptions

The implementation follows the main relative-quantification sequence
described by Taylor *et al.*:

| Procedure in Taylor *et al.*                                                                       | Current implementation                                               |
|----------------------------------------------------------------------------------------------------|----------------------------------------------------------------------|
| Average technical-replicate Cq values for each biological sample and target                        | Implemented for repeated `Sample`-`Target` pairs                     |
| Calculate relative quantity from the control-group mean Cq and each sample Cq                      | Implemented                                                          |
| Use assay-specific PCR efficiency in the relative-quantity calculation                             | Currently assumes 100% efficiency and uses base `2` for every target |
| Normalize each sample using the geometric mean of two or more reference-target relative quantities | Implemented; users select reference genes with `controls`            |
| Calculate group-level normalized expression                                                        | Implemented as the geometric mean                                    |
| Perform inference on log-transformed normalized expression per biological sample                   | Implemented by default with `stats_value = "log2_norm_exp"`          |

The package warns about technical-replicate variability but does not
currently perform inter-run calibration, amplification-efficiency
estimation, automatic outlier removal, or reference-gene stability
validation. Therefore, results correspond to the paper’s procedure only
when the remaining experimental and preprocessing requirements have
already been addressed and the assumption of approximately 100%
amplification efficiency is appropriate.

Save the analyzed table by supplying both arguments:

``` r
analysis <- analyze_qPCR(
  path_to_qPCR_res = "data/quantstudio_results.csv",
  path_to_group_label = "data/sample_metadata.csv",
  out_dir_path = "results",
  save_results = TRUE
)
```

### Run statistics

``` r
statistical_result <- stats_qPCR(
  analysis,
  stats_type = "ttest",
  stats_value = "log2_norm_exp"
)
```

Available expression columns are `"log2_norm_exp"` and `"norm_exp"`.
`"log2_norm_exp"` is the default and is strongly recommended for the
package’s parametric t-tests and ANOVA. Relative-expression ratios are
positive, multiplicative measurements and are often right-skewed or
approximately log-normally distributed on the linear scale. Log2
transformation makes fold changes additive and generally better supports
the normality and constant-variance assumptions of these tests.

Log transformation does not guarantee that the statistical assumptions
are met. When possible, inspect the model residuals and consider a more
appropriate statistical method if substantial violations remain. Use
`"norm_exp"` for inference only when its use is justified for the data
and selected model. Untestable targets, such as targets with fewer than
two observations in a group, are skipped with a warning.

### Select plot comparisons

The `comparison` argument controls which calculated comparisons are
displayed:

``` r
plot_qPCR(
  statistical_result,
  plot_individual = TRUE,
  comparison = "control"
)
```

Available modes are:

- `"all"`: display every calculated comparison;
- `"control"`: display comparisons involving the reference group;
- `"significant"`: display comparisons passing the significance
  threshold;
- `"custom"`: display pairs supplied through `stats_group`.

Each custom entry must contain exactly two group names:

``` r
plot_qPCR(
  statistical_result,
  plot_individual = TRUE,
  comparison = "custom",
  stats_group = list(
    c("Control", "Treatment_A"),
    c("Treatment_A", "Treatment_B")
  )
)
```

Group names must exactly match the `Group` column.

## Plot customization

### Named colors

Named color vectors are matched by group name, so their order does not
matter:

``` r
group_colors <- c(
  Treatment_B = "#4DBBD5",
  Control = "#FFFFFF",
  Treatment_A = "#E64B35"
)

plot_qPCR(
  analysis,
  palette = group_colors,
  plot_individual = TRUE
)
```

Unnamed color vectors are assigned positionally according to
`Group_order`. Built-in palette names include `"grouped"`, `"viridis"`,
`"npg"`, `"aaas"`, `"lancet"`, `"jama"`, `"nejm"`, `"jco"`, and several
Khroma palettes.

### Reproducible jitter

Individual plots jitter sample points horizontally by default. The
default seed is `1`, so repeated plotting produces the same point
positions:

``` r
plot_qPCR(
  analysis,
  plot_individual = TRUE,
  jitter = TRUE,
  jitter_width = 0.1,
  jitter_seed = 1
)
```

Disable jitter to center the points over each group:

``` r
plot_qPCR(
  analysis,
  plot_individual = TRUE,
  jitter = FALSE
)
```

## Saving results

The complete workflow can save:

- `analyzed_qPCR_res.tsv`;
- `qPCR_stats.tsv`;
- `qPCR_plot.tiff` and `qPCR_plot.svg`;
- one TIFF and SVG file per target when `plot_individual = TRUE`.

Control the output stages independently:

``` r
result <- qPCR(
  path_to_qPCR_res = "data/quantstudio_results.csv",
  path_to_group_label = "data/sample_metadata.csv",
  out_dir_path = "results",
  stats_type = "ttest",
  save_analysis = TRUE,
  save_stats = TRUE,
  save_plots = TRUE
)
```

When calling `plot_qPCR()` directly, plots are returned without being
saved by default. Use `save_tiff = TRUE` and/or `save_svg = TRUE`
together with `out_dir_path` to write files.

## Reference

Taylor SC, Nadeau K, Abbasi M, Lachance C, Nguyen M, Fenrich J. The
Ultimate qPCR Experiment: Producing Publication Quality, Reproducible
Data the First Time. *Trends in Biotechnology*. 2019;37(7):761–774.
[doi:10.1016/j.tibtech.2018.12.002](https://doi.org/10.1016/j.tibtech.2018.12.002).
