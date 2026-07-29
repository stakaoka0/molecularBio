# classes.R

#' Print qpcr analysis
#'
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
