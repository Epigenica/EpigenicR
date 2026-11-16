#' Stacked bar plot of feature-set composition per sample, in percent
#'
#' For each marker and sample, plots how the signal is distributed over the
#' feature sets (experiments in \code{epk$mse}) as a percentage. Each bar adds
#' up to 100\%. One panel is drawn per marker.
#'
#' @param epk An \code{EPK} object.
#' @param markers Character vector of assay names (e.g. \code{"H3K27ac"}).
#'   Default \code{NULL} uses the markers present in every selected feature
#'   set.
#' @param features Character vector of experiment names in \code{epk$mse}.
#'   Default \code{NULL} uses all experiments. For a genomic distribution,
#'   pass the sets from \code{\link{build_genomic_feature_sets}}, e.g.
#'   \code{c("tss", "promoter", "utr5", "utr3", "exon", "intron",
#'   "intergenic")}.
#' @param samples Character vector of sample (column) names to keep.
#'   Default \code{NULL} keeps all.
#' @param metric What the percentage is based on:
#'   \describe{
#'     \item{\code{"loci"}}{(default) the number of loci in each feature set
#'       with signal above \code{threshold}.}
#'     \item{\code{"signal"}}{the total signal in each feature set (mean
#'       signal of each locus multiplied by its width, summed). With
#'       non-overlapping feature sets that cover the genome, this is the
#'       genomic distribution of the signal.}
#'   }
#' @param threshold Numeric; with \code{metric = "loci"}, a locus counts as
#'   marked if its value is greater than this. Default \code{1} (above the
#'   genome-wide average for RPGC). Ignored for \code{metric = "signal"}.
#'@param show_input Logical; add an \code{"Input (median)"} bar to each panel.
#'   The distribution over the feature sets is calculated for every INPUT
#'   sample in the same way as for the markers, and the bar shows the median
#'   over the INPUT samples. Needs the INPUT marker in the EPK. Default
#'   \code{FALSE}.
#' @param input_marker Name of the INPUT assay. Default \code{"INPUT"}.
#' @param ncol Number of panel columns. Default \code{NULL} puts all markers
#'   in one row.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The percentage is each feature set's value divided by the sum over all
#' selected feature sets, per sample and marker. With overlapping feature
#' sets of different sizes (e.g. gene bodies, CpG islands) this shows the
#' composition of marked loci, not a genomic distribution. For a genomic
#' distribution, use the non-overlapping sets from
#' \code{\link{build_genomic_feature_sets}} together with
#' \code{metric = "signal"}.
#'
#' @examples
#' \dontrun{
#' plot_features_stacked_bar(epk)
#' plot_features_stacked_bar(epk, markers = "H3K27ac", threshold = 2)
#'
#' # Genomic distribution of signal over TSS, UTRs, exons, introns, intergenic
#' sets <- build_genomic_feature_sets("data/gencode.v38.annotation.gtf",
#'                                    chrom_sizes = "path/to/a.bw")
#' epk <- add_features_to_epk(epk, annotations = sets,
#'                            pipeline_output_path = "/path/to/project")
#' plot_features_stacked_bar(epk, features = names(sets), metric = "signal",
#'                           show_input = TRUE)
#' }
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap scale_y_continuous
#'   theme_classic labs theme element_text
#' @export
plot_features_stacked_bar <- function(epk, markers = NULL, features = NULL,
                                      samples = NULL,
                                      metric = c("loci", "signal"),
                                      threshold = 1,
                                      show_input = FALSE,
                                      input_marker = "INPUT",
                                      ncol = NULL) {
  metric <- match.arg(metric)

  if (is.null(epk) || !is.list(epk) || is.null(epk$mse)) {
    stop("'epk' must be a valid EPK object containing an 'mse' slot.")
  }

  all_features <- names(epk$mse)
  if (is.null(features)) features <- all_features
  unknown <- setdiff(features, all_features)
  if (length(unknown) > 0) {
    stop("'features' not found in epk$mse: ", paste(unknown, collapse = ", "))
  }

  assay_names <- lapply(features, function(ft) {
    SummarizedExperiment::assayNames(epk$mse[[ft]])
  })
  common_markers <- Reduce(intersect, assay_names)
  if (is.null(markers)) markers <- setdiff(common_markers, input_marker)
  missing_markers <- setdiff(markers, common_markers)
  if (length(missing_markers) > 0) {
    stop(
      "Marker(s) not present in every selected feature set: ",
      paste(missing_markers, collapse = ", ")
    )
  }

  # Locus widths are needed to turn mean signal into total signal
  widths <- NULL
  if (metric == "signal") {
    widths <- lapply(features, function(ft) {
      se <- epk$mse[[ft]]
      if (!methods::is(se, "RangedSummarizedExperiment")) {
        stop(
          "metric = 'signal' needs genomic ranges, but experiment '", ft,
          "' has none."
        )
      }
      as.numeric(GenomicRanges::width(SummarizedExperiment::rowRanges(se)))
    })
    names(widths) <- features
  }

  plot_df <- do.call(rbind, lapply(markers, function(mk) {
    do.call(rbind, lapply(features, function(ft) {
      mat <- SummarizedExperiment::assay(epk$mse[[ft]], mk)
      value <- if (metric == "signal") {
        colSums(mat * widths[[ft]], na.rm = TRUE)
      } else {
        colSums(mat > threshold, na.rm = TRUE)
      }
      data.frame(
        marker  = mk,
        sample  = colnames(mat),
        feature = ft,
        value   = value,
        row.names = NULL
      )
    }))
  }))

  if (!is.null(samples)) {
    unknown_samples <- setdiff(samples, plot_df$sample)
    if (length(unknown_samples) > 0) {
      stop("'samples' not found: ", paste(unknown_samples, collapse = ", "))
    }
    plot_df <- plot_df[plot_df$sample %in% samples, , drop = FALSE]
  }
  sample_levels <- unique(plot_df$sample)

  # Reference bar: median distribution of the INPUT samples
  if (isTRUE(show_input)) {
    if (!input_marker %in% common_markers) {
      stop(
        "show_input = TRUE needs the '", input_marker, "' marker in every ",
        "selected feature set. Build the EPK with markers_to_exclude = character(0)."
      )
    }
    input_label <- "Input (median)"
    input_df <- do.call(rbind, lapply(features, function(ft) {
      mat <- SummarizedExperiment::assay(epk$mse[[ft]], input_marker)
      value <- if (metric == "signal") {
        colSums(mat * widths[[ft]], na.rm = TRUE)
      } else {
        colSums(mat > threshold, na.rm = TRUE)
      }
      data.frame(sample = colnames(mat), feature = ft, value = value,
                 row.names = NULL)
    }))
    if (!is.null(samples)) {
      input_df <- input_df[input_df$sample %in% samples, , drop = FALSE]
    }
    # Percentage per INPUT sample, then the median over the INPUT samples
    input_total <- stats::ave(input_df$value, input_df$sample, FUN = sum)
    input_df <- input_df[input_total > 0, , drop = FALSE]
    input_total <- input_total[input_total > 0]
    if (nrow(input_df) == 0) {
      stop("No '", input_marker, "' data found for the selected samples.")
    }
    input_df$percent <- 100 * input_df$value / input_total
    input_median <- tapply(input_df$percent, input_df$feature, stats::median)
    input_bar <- do.call(rbind, lapply(markers, function(mk) {
      data.frame(
        marker  = mk,
        sample  = input_label,
        feature = features,
        value   = as.numeric(input_median[features]),
        row.names = NULL
      )
    }))
    plot_df <- rbind(plot_df, input_bar)
    sample_levels <- c(sample_levels, input_label)
  }

  # Percentage per sample and marker
  totals <- stats::ave(plot_df$value, plot_df$marker, plot_df$sample, FUN = sum)
  plot_df$percent <- ifelse(totals > 0, 100 * plot_df$value / totals, 0)

  plot_df$feature <- factor(plot_df$feature, levels = features)
  plot_df$marker  <- factor(plot_df$marker, levels = markers)
  plot_df$sample  <- factor(plot_df$sample, levels = sample_levels)
  if (is.null(ncol)) ncol <- length(markers)

  y_label <- if (metric == "signal") "Share of signal (%)" else "Frequency (%)"
  plot_title <- if (metric == "signal") {
    "Distribution of signal across feature sets"
  } else {
    paste0("Loci with signal > ", threshold, " per feature set")
  }

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$sample, y = .data$percent, fill = .data$feature)
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~marker, ncol = ncol) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100.001), expand = c(0, 0),
      breaks = seq(0, 100, 20)
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Sample",
      y = y_label,
      fill = "Feature set",
      title = plot_title
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
