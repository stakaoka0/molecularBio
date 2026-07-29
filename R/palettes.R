# resolve_palette <- function(
#   palette,
#   groups
# ) {
#   n <- length(groups)
#   cols <- switch(
#     palette,
#     nature =
#       ggsci::pal_npg()(
#         max(n, 10)
#       )[seq_len(n)],
#     nejm =
#       ggsci::pal_nejm()(
#         max(n, 8)
#       )[seq_len(n)],
#     lancet =
#       ggsci::pal_lancet()(
#         max(n, 9)
#       )[seq_len(n)],
#     jama =
#       ggsci::pal_jama()(
#         max(n, 7)
#       )[seq_len(n)],
#     viridis =
#       viridisLite::viridis(n),
#     stop("Unknown palette")
#   )
#   stats::setNames(
#     cols,
#     groups
#   )
# }
