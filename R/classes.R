# classes.R

#' Print qpcr analysis
#'
#' @param x A `qpcr_analysis` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return `x`, invisibly.
#' @export

print.qpcr_analysis <- function(
    x,
    ...
) {
  cat("\n")
  cat("qPCR Analysis Object\n")
  cat("--------------------\n")
  cat(
    "Reference group:",
    x$reference_group,
    "\n"
  )
  cat(
    "Control genes:",
    paste(
      x$controls,
      collapse = ", "
    ),
    "\n"
  )
  cat(
    "Targets:",
    length(
      unique(
        x$data$Target
      )
    ),
    "\n"
  )
  invisible(x)
}


#' Summary qpcr analysis
#'
#' @param object A `qpcr_analysis` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A one-row tibble summarizing the analysis.
#' @export

summary.qpcr_analysis <- function(
    object,
    ...
) {
  tibble::tibble(
    targets =
      length(
        unique(
          object$data$Target
        )
      ),
    groups =
      length(
        unique(
          object$data$Group
        )
      ),
    controls =
      paste(
        object$controls,
        collapse = ", "
      ),
    reference_group =
      object$reference_group
  )
}


#' Plot method
#'
#' @param x A `qpcr_analysis` object.
#' @param ... Arguments passed to [plot_qPCR()].
#'
#' @return Invisibly, the plots returned by [plot_qPCR()].
#' @export

plot.qpcr_analysis <- function(
    x,
    ...
) {
  plot_qPCR(
    object = x,
    ...
  )
}


#' Print statistical results
#'
#' @param x A `qpcr_stats` object.
#' @param ... Additional arguments passed to the next print method.
#'
#' @return `x`, invisibly.
#' @export

print.qpcr_stats <- function(
    x,
    ...
) {
  NextMethod()
  cat(
    "\nStatistical results available:\n"
  )
  cat(
    nrow(x$stats),
    "comparisons\n"
  )
  invisible(x)
}
