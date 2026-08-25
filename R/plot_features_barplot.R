#' Bar plot of loci with signal per feature set and marker
#'
#' For each experiment (feature/annotation set) in an \code{EPK} object,
#' counts how many loci have non-zero total signal for each marker, and
#' plots the result as a stacked bar chart.
#'
#' @param epk An \code{EPK} object.
#' @param features Character vector of experiment names in \code{epk$mse} to
#'   include. Default \code{NULL} uses every experiment
#'   (\code{names(epk$mse)}).
#' @param markers Character vector of marker/assay names to include. Default
#'   \code{NULL} uses every assay present in the first selected experiment.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' For each \code{(feature, marker)} pair, a locus counts as "with signal" if
#' its total (row-summed) value across samples is greater than zero
#' (\code{NA} values are treated as \code{0} via \code{na.rm = TRUE}, so an
#' all-\code{NA} locus is indistinguishable from a genuinely zero-signal
#' locus). Bars use \code{position = "fill"}, so each marker's bar shows
#' every included feature set's share of the \emph{combined} loci-with-signal
#' count for that marker across feature sets \eqn{-} not each feature set's
#' own percentage of its own total loci. Feature sets with very different
#' total loci counts are therefore not directly comparable to each other by
#' segment height alone.
#'
#' @examples
#' \dontrun{
#' plot_features_barplot(epk)
#'
#' plot_features_barplot(
#'   epk,
#'   features = c("genes", "cpgs"),
#'   markers  = c("H3K4me3", "H3K27ac")
#' )
#' }
#'
#' @importFrom SummarizedExperiment assays assayNames
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_bar theme_classic labs theme element_text
#' @export
plot_features_barplot <- function(epk, features = NULL, markers = NULL) {
  if (is.null(epk) || !is.list(epk) || is.null(epk$mse)) {
    stop("'epk' must be a valid EPK object containing an 'mse' slot.")
  }

  all_features <- names(epk$mse)
  if (is.null(features)) {
    features <- all_features
  } else {
    unknown_features <- setdiff(features, all_features)
    if (length(unknown_features) > 0) {
      stop(
        "'features' contains name(s) not found in epk$mse: ",
        paste(unknown_features, collapse = ", ")
      )
    }
  }

  if (is.null(markers)) {
    markers <- SummarizedExperiment::assayNames(epk$mse[[features[1]]])
  }

  feature_freq <- setNames(
    data.frame(
      matrix(NA_real_, nrow = length(features), ncol = length(markers)),
      row.names = features
    ),
    markers
  )

  for (ft in features) {
    ft_assay_names <- SummarizedExperiment::assayNames(epk$mse[[ft]])
    missing_markers <- setdiff(markers, ft_assay_names)
    if (length(missing_markers) > 0) {
      stop(
        "Experiment '", ft, "' is missing marker(s): ",
        paste(missing_markers, collapse = ", ")
      )
    }

    ft_assays <- SummarizedExperiment::assays(epk$mse[[ft]])[markers]
    feature_freq[ft, markers] <- sapply(
      ft_assays,
      function(x) sum(rowSums(x, na.rm = TRUE) > 0, na.rm = TRUE)
    )
  }

  plot_df <- feature_freq |>
    tibble::rownames_to_column("feature") |>
    tidyr::pivot_longer(-"feature", names_to = "marker", values_to = "count")

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$marker, y = .data$count, fill = .data$feature)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::theme_classic() +
    ggplot2::labs(y = "Percentage of features with signal") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
