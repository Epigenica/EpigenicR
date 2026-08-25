#' Plot fragment (library complexity) distribution QC across samples
#'
#' For each marker, ranks samples by a chosen library-complexity metric and
#' flags samples whose value falls well below that marker's own median,
#' highlighting outliers on a faceted rank plot.
#'
#' @param data Backward-compatible alias for \code{stats_summary}.
#' @param epk Optional EPK object. If supplied, \code{epk$tables$stats_summary}
#'   is used as plotting input.
#' @param stats_summary Optional data frame containing QC statistics. If
#'   provided, it takes precedence over \code{epk} and \code{data}. Must
#'   include \code{marker}, \code{unique_fragments_col}, and
#'   \code{sample_labeling}.
#' @param markers Character vector of markers to include. Default \code{NULL}
#'   uses every marker present in the data.
#' @param unique_fragments_col Character; name of the column holding the
#'   library-complexity metric to rank samples on (e.g. \code{"library_size"}
#'   or \code{"final_mapped"} \eqn{-} pick whichever your pipeline reports as
#'   unique fragments). Default \code{"library_size"}.
#' @param sample_labeling Character; column used to label flagged samples.
#'   One of \code{"sample_id_rep"}, \code{"sample_id"}, or \code{"map_id"}.
#'   Default \code{"sample_id_rep"}.
#' @param warning_threshold Numeric in (0, 1); a sample is flagged
#'   \code{"warning"} when its \code{unique_fragments_col} value falls below
#'   this fraction of its marker's median. Default \code{0.50}.
#' @param critical_threshold Numeric in (0, \code{warning_threshold}); samples
#'   below this fraction of the marker median are flagged \code{"critical"}.
#'   Default \code{0.25}.
#' @param marker_levels Optional character vector defining facet order. If
#'   \code{NULL} (default), uses the order found in the data.
#' @param ncol Integer; number of facet columns. Default \code{2}.
#' @param label_flagged Logical; if \code{TRUE} (default), non-\code{"good"}
#'   samples are labeled with \code{sample_labeling}. Uses \pkg{ggrepel} when
#'   installed to avoid overlapping labels, otherwise falls back to
#'   \code{ggplot2::geom_text()}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{plot}: a \code{ggplot} object, faceted by marker.
#'   \item \code{data}: the per-sample data frame used to build the plot,
#'     with added \code{marker_median}, \code{ratio_to_median},
#'     \code{qc_class} (\code{"good"}/\code{"warning"}/\code{"critical"}), and
#'     \code{rank} columns.
#' }
#'
#' @details
#' Thresholds are relative to each marker's own median, not a fixed
#' vendor pass/fail cutoff \eqn{-} a marker with generally lower yield does
#' not systematically flag more samples just because of its baseline depth.
#' Rows with \code{NA} in \code{unique_fragments_col} get \code{qc_class = NA}
#' rather than being swept into \code{"critical"}.
#'
#' @examples
#' \dontrun{
#' qc <- plot_fragment_distribution(epk = epk)
#' qc$plot
#' qc$data
#'
#' qc <- plot_fragment_distribution(
#'   epk = epk,
#'   markers = c("H3K4me3", "H3K27ac", "H3K27me3", "5mC"),
#'   unique_fragments_col = "final_mapped"
#' )
#' }
#'
#' @importFrom dplyr filter group_by mutate ungroup arrange desc distinct case_when row_number
#' @importFrom ggplot2 ggplot aes geom_hline geom_point geom_text facet_wrap
#' @importFrom ggplot2 scale_y_log10 scale_color_manual scale_size_manual
#' @importFrom ggplot2 labs theme_minimal theme element_blank element_text
#' @importFrom stats median
#' @export
plot_fragment_distribution <- function(
  data = NULL,
  epk = NULL,
  stats_summary = NULL,
  markers = NULL,
  unique_fragments_col = "library_size",
  sample_labeling = c("sample_id_rep", "sample_id", "map_id"),
  warning_threshold = 0.50,
  critical_threshold = 0.25,
  marker_levels = NULL,
  ncol = 2,
  label_flagged = TRUE
) {
  sample_labeling <- match.arg(sample_labeling)

  if (!(critical_threshold > 0 && critical_threshold < warning_threshold && warning_threshold < 1)) {
    stop("Require 0 < critical_threshold < warning_threshold < 1.")
  }

  data <- .resolve_stats_summary_input(
    data = data,
    epk = epk,
    stats_summary = stats_summary,
    fn_name = "plot_fragment_distribution"
  )

  required_cols <- c("marker", unique_fragments_col, sample_labeling)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ", paste(missing_cols, collapse = ", "),
      ". Available numeric columns: ",
      paste(names(data)[vapply(data, is.numeric, logical(1))], collapse = ", ")
    )
  }

  if (!is.null(markers)) {
    unknown_markers <- setdiff(markers, unique(data$marker))
    if (length(unknown_markers) > 0) {
      stop(
        "'markers' contains value(s) not found in data$marker: ",
        paste(unknown_markers, collapse = ", ")
      )
    }
    data <- dplyr::filter(data, .data$marker %in% markers)
  }

  if (nrow(data) == 0) {
    stop("No rows remain after filtering to 'markers'.")
  }

  if (is.null(marker_levels)) {
    marker_levels <- unique(data$marker)
  }
  data$marker <- factor(data$marker, levels = marker_levels)

  qc_data <- data |>
    dplyr::group_by(.data$marker) |>
    dplyr::mutate(
      marker_median   = stats::median(.data[[unique_fragments_col]], na.rm = TRUE),
      ratio_to_median = .data[[unique_fragments_col]] / .data$marker_median,
      qc_class = dplyr::case_when(
        is.na(.data$ratio_to_median)                 ~ NA_character_,
        .data$ratio_to_median >= warning_threshold    ~ "good",
        .data$ratio_to_median >= critical_threshold   ~ "warning",
        TRUE                                          ~ "critical"
      )
    ) |>
    dplyr::arrange(.data$marker, dplyr::desc(.data[[unique_fragments_col]])) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup()

  ref_lines <- qc_data |>
    dplyr::distinct(.data$marker, .data$marker_median) |>
    dplyr::mutate(
      warn_line = .data$marker_median * warning_threshold,
      crit_line = .data$marker_median * critical_threshold
    )

  qc_colors <- c(good = "#2a78d6", warning = "#fab219", critical = "#d03b3b")
  qc_data$qc_class <- factor(qc_data$qc_class, levels = names(qc_colors))

  p <- ggplot2::ggplot(qc_data, ggplot2::aes(x = .data$rank, y = .data[[unique_fragments_col]])) +
    ggplot2::geom_hline(
      data = ref_lines, ggplot2::aes(yintercept = .data$marker_median),
      color = "grey40", linewidth = 0.6, alpha = 0.6
    ) +
    ggplot2::geom_hline(
      data = ref_lines, ggplot2::aes(yintercept = .data$warn_line),
      color = "#fab219", linewidth = 0.6, linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      data = ref_lines, ggplot2::aes(yintercept = .data$crit_line),
      color = "#d03b3b", linewidth = 0.6, linetype = "dashed"
    ) +
    ggplot2::geom_point(ggplot2::aes(color = .data$qc_class, size = .data$qc_class), alpha = 0.75) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(values = qc_colors, name = "QC flag", na.translate = FALSE) +
    ggplot2::scale_size_manual(
      values = c(good = 1.4, warning = 2, critical = 2.4), guide = "none", na.value = 1.4
    ) +
    ggplot2::facet_wrap(~marker, scales = "free_x", ncol = ncol) +
    ggplot2::labs(
      x = "Libraries ranked by unique fragments (high to low)",
      y = paste0(unique_fragments_col, " (log10)"),
      title = "Library complexity QC — unique fragments per sample",
      subtitle = "Thresholds are relative to each marker's own median"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (isTRUE(label_flagged)) {
    flagged <- dplyr::filter(qc_data, !is.na(.data$qc_class) & .data$qc_class != "good")
    if (nrow(flagged) > 0) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = flagged,
          ggplot2::aes(label = .data[[sample_labeling]], color = .data$qc_class),
          size = 2.6, max.overlaps = Inf, show.legend = FALSE
        )
      } else {
        p <- p + ggplot2::geom_text(
          data = flagged,
          ggplot2::aes(label = .data[[sample_labeling]], color = .data$qc_class),
          size = 2.6, vjust = -0.6, show.legend = FALSE
        )
      }
    }
  }

  list(plot = p, data = qc_data)
}
