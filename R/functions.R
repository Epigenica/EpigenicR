#' Plot beta-value density distributions from EPK
#'
#' Plots density distributions of beta values (raw, quantile-normalized, or both)
#' for a given feature in an EPK object. Facets by normalization method if both.
#'
#' @param epk An EPK object.
#' @param feature Character; experiment name (e.g. "cpg_islands", "union_peaks").
#' @param method Character; one of "raw", "quantile", or "both" (default: "both").
#' @param assay_raw Character; name of the raw beta assay (default: "beta_raw").
#' @param assay_qn Character; name of the quantile-normalized beta assay (default: "beta_qn").
#' @return A ggplot object (density plot, faceted if both).
## Example usage:
#' @examples
#' # Assuming you have an EPK object named epk
#' plot_beta_density(epk, feature = "cpg_islands", method = "both")
#' plot_beta_density(epk, feature = "genes", method = "raw")
#' @export
plot_beta_density <- function(epk,
                              feature = "cpg_islands",
                              method = c("both", "raw", "quantile"),
                              assay_raw = "beta_raw",
                              assay_qn = "beta_qn") {
  method <- match.arg(method)
  se <- MultiAssayExperiment::experiments(epk$mse)[[feature]]
  dfs <- list()
  if (method %in% c("raw", "both")) {
    if (!assay_raw %in% names(SummarizedExperiment::assays(se))) {
      stop(sprintf("Assay '%s' not found in experiment '%s'", assay_raw, feature))
    }
    df_raw <- as.data.frame(SummarizedExperiment::assay(se, assay_raw)) %>%
      tibble::rownames_to_column("region") %>%
      tidyr::pivot_longer(-region, names_to = "sample", values_to = "beta") %>%
      dplyr::mutate(method = "Before QN")
    dfs[["raw"]] <- df_raw
  }
  if (method %in% c("quantile", "both")) {
    if (!assay_qn %in% names(SummaraturExperiment::assays(se))) {
      stop(sprintf("Assay '%s' not found in experiment '%s'", assay_qn, feature))
    }
    df_qn <- as.data.frame(SummarizedExperiment::assay(se, assay_qn)) %>%
      tibble::rownames_to_column("region") %>%
      tidyr::pivot_longer(-region, names_to = "sample", values_to = "beta") %>%
      dplyr::mutate(method = "After QN")
    dfs[["qn"]] <- df_qn
  }
  df_plot <- dplyr::bind_rows(dfs)
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = beta, colour = sample)) +
    ggplot2::geom_density(linewidth = 0.7, alpha = 0.8) +
    ggplot2::facet_wrap(~method, ncol = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Beta-value distributions before and after quantile normalization",
      x = "Beta value",
      y = "Density"
    )
  p
}
#' Extract marker names from sample IDs
#' @param id A character string representing the sample ID.
#' @param markers A vector of marker names to search for within the sample ID.
#' @return A character vector containing the extracted marker names.
## Example usage:
#' @examples
#' # Compute and add sample correlations for all experiments
#' epk <- compute_sample_cor(epk, method = "pearson", transform = "none")
#' # Access correlations:
#' # epk$derived$sample_cor
#' @export
#' @examples
#' id <- "Sample1_H3K4me3_Rep1"
#' markers <- c("INPUT","H3K4me3","H3K9me3","H3K9ac","5mC","CXXC")
#' extract_marker_names(id, markers)
#'
#' @importFrom stringr str_extract str_c
extract_marker_names <- function(id, markers){
  if ( is.integer(id)){
    print("The ID should be a character string.")
  }
  extrc_markers <- stringr::str_match(
    id,
    "^([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_(rep[0-9]+|pooled)\\.(hg38)\\.(scaled|unscaled)\\.bw$")[4]
  return(extrc_markers)
}

#' Resolve QC input to a stats_summary-like data frame
#' @keywords internal
.resolve_stats_summary_input <- function(data = NULL, epk = NULL, stats_summary = NULL, fn_name = "plot_qc_stats") {
  provided <- c(
    data = !is.null(data),
    epk = !is.null(epk),
    stats_summary = !is.null(stats_summary)
  )

  if (!any(provided)) {
    stop("Provide one of 'stats_summary', 'epk', or 'data' to ", fn_name, "().")
  }

  if (sum(provided) > 1) {
    warning(
      "Multiple input sources provided (",
      paste(names(provided)[provided], collapse = ", "),
      "). Using priority: stats_summary > data > epk."
    )
  }

  if (!is.null(stats_summary)) {
    data_in <- stats_summary
  } else if (!is.null(data)) {
    data_in <- data
  } else {
    if (!is.list(epk) || is.null(epk$tables) || is.null(epk$tables$stats_summary)) {
      stop("'epk' must contain epk$tables$stats_summary.")
    }
    data_in <- epk$tables$stats_summary
  }

  if (!is.data.frame(data_in)) {
    data_in <- tryCatch(
      as.data.frame(data_in, stringsAsFactors = FALSE),
      error = function(e) {
        stop("Resolved stats summary input is not a data frame and cannot be coerced.")
      }
    )
  }

  data_in
}

#' Plot QC statistics across markers, conditions, and replicates
#'
#' Generates boxplots with overlaid jittered points for QC metrics across markers
#' and experimental conditions. Can return either a static multi-panel figure
#' composed with \pkg{patchwork} or an interactive \pkg{plotly} subplot grid.
#'
#' @param data Backward-compatible alias for \code{stats_summary}. A data frame
#'   containing QC statistics.
#'
#' @param epk Optional EPK object. If supplied,\code{epk$tables$stats_summary}
#'   is used as plotting input.
#'
#' @param stats_summary Optional data frame containing QC statistics. If
#'   provided, it takes precedence over \code{epk} and \code{data}.
#'   Must include a \code{marker} column and (for full functionality)
#'   \code{replicate}. If \code{condition} is missing, all rows are treated
#'   as a single condition. All numeric columns are treated as candidate QC
#'   statistics unless \code{stats} is provided.
#'
#' @param condition Character vector specifying which conditions to plot.
#'   \itemize{
#'     \item \code{NULL}: plot all available conditions as-is (default)
#'     \item \code{"All"}: collapse all rows into a single condition named \code{"All"}
#'     \item otherwise: one or more condition names to subset
#'   }
#'
#' @param stats Character vector of numeric column names to plot. If \code{NULL}
#'   (default), all numeric columns in \code{data} are used.
#'
#' @param stats_exclude Optional character vector of numeric column names to
#'   exclude from plotting. Applied after \code{stats} selection. Useful when
#'   \code{stats = NULL} but a few numeric columns should be omitted.
#'
#' @param marker_levels Optional character vector defining the order of markers
#'   on the x-axis. If \code{NULL} (default), uses the order found in
#'   \code{data$marker}.
#'
#' @param legend_position Character string indicating where to place the legend.
#'   For \code{engine = "ggplot"}, this is passed to
#'   \code{theme(legend.position = ...)} (e.g. \code{"bottom"}, \code{"right"}).
#'   For \code{engine = "plotly"}, \code{"bottom"} produces a horizontal legend,
#'   otherwise a vertical legend is used.
#'
#' @param save_plots Logical; if \code{TRUE} and \code{engine = "ggplot"}, saves
#'   each composed condition-level figure as a PNG file. Saving interactive
#'   plotly output is not handled by this function.
#'
#' @param save_dir Character path to directory where plots will be saved.
#'   Required only if \code{save_plots = TRUE} and \code{engine = "ggplot"}.
#'
#' @param ncol Integer; number of columns used when composing plots per condition.
#'   Used to determine which panels are on the bottom row (to show x-axis labels
#'   only there) and should match the layout used by \code{patchwork::wrap_plots()}
#'   or the \code{plotly::subplot()} grid.
#'
#' @param sample_labeling Character; column in \code{data} used to color points.
#'   One of \code{"map_id"}, \code{"sample_id_rep"}, \code{"sample_id"},
#'   or \code{"replicate"}. Default: \code{"map_id"}.
#'
#' @param engine Character; output engine. One of \code{"ggplot"} or \code{"plotly"}.
#'   \itemize{
#'     \item \code{"ggplot"}: returns a static multi-panel patchwork/ggplot object
#'       with collected guides (single legend) and x-axis labels shown only on the
#'       bottom row.
#'     \item \code{"plotly"}: returns an interactive plotly subplot grid created
#'       via \code{plotly::ggplotly()} + \code{plotly::subplot()}. The legend is
#'       displayed only once (first panel) to avoid duplicated entries.
#'   }
#'
#' @return If \code{condition} resolves to a single condition, returns a single
#'   object:
#'   \itemize{
#'     \item for \code{engine = "ggplot"}: a patchwork/ggplot object
#'     \item for \code{engine = "plotly"}: a plotly htmlwidget
#'   }
#'   If multiple conditions are requested, returns a named list of such objects,
#'   one per condition (names correspond to condition values).
#'
#' @details
#' Each panel consists of:
#' \itemize{
#'   \item boxplots grouped by marker (outliers hidden),
#'   \item jittered points showing individual replicates,
#'   \item fill mapped to marker and color mapped to replicate.
#' }
#' For multi-panel figures, x-axis tick labels are hidden for all panels except
#' those in the bottom row (computed from \code{ncol}).
#' The \code{condition} argument is primarily a subsetting selector; use
#' \code{"All"} only when you explicitly want to collapse across conditions.
#'
#' @examples
#' \dontrun{
#' # Static (publication-style) figure with one legend at the bottom
#' p <- plot_qc_stats(
#'   data = qc_df,
#'   condition = "Ctrl",
#'   stats = c("final_mapped", "library_size"),
#'   marker_levels = c("INPUT", "H3K4me3", "H3K9me3", "H3K9ac", "5mC", "CXXC"),
#'   engine = "ggplot",
#'   legend_position = "bottom",
#'   ncol = 3
#' )
#' p
#'
#' # Interactive plotly grid (legend shown once)
#' fig <- plot_qc_stats(
#'   data = qc_df,
#'   condition = "Ctrl",
#'   engine = "plotly",
#'   legend_position = "bottom",
#'   sample_labeling = "map_id",
#'   ncol = 3
#' )
#' fig
#'
#' # EPK-centered QC plotting
#' qc_plots <- plot_qc_stats(
#'   epk = epk,
#'   legend_position = "none",
#'   save_plots = TRUE,
#'   save_dir = results_qc,
#'   ncol = 3,
#'   engine = "plotly",
#'   sample_labeling = "map_id"
#' )
#'
#' # Direct stats_summary input
#' qc_plots_tbl <- plot_qc_stats(
#'   stats_summary = qc_df,
#'   engine = "ggplot",
#'   ncol = 3
#' )
#'
#' # Batch-aware workflow (derive metadata outside plot_qc_stats)
#' stats_summary <- stats_summary %>%
#'   dplyr::mutate(
#'     create_metadata_df(map_id_vector = stats_summary$map_id, bw_files = NULL)[,
#'       c("marker", "batch", "sample_id", "replicate")
#'     ],
#'     sample_id_rep = paste0(sample_id, "_", replicate),
#'     sample_id_rep_batch = paste0(sample_id_rep, "_", batch)
#'   )
#'
#' qc_plots_A1 <- plot_qc_stats(
#'   stats_summary %>% dplyr::select(-msr) %>% dplyr::filter(batch == "A1"),
#'   legend_position = "none",
#'   save_plots = TRUE,
#'   save_dir = results_qc,
#'   ncol = 3,
#'   engine = "plotly",
#'   sample_labeling = "sample_id_rep_batch"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot geom_boxplot geom_point theme_bw labs guides theme
#' @importFrom ggplot2 element_text element_blank element_line position_jitter ggsave
#' @importFrom dplyr filter select
#' @importFrom tidyselect where
#' @export
plot_qc_stats <- function(
  data = NULL,
  epk = NULL,
  stats_summary = NULL,
    condition = NULL,
    stats = NULL,
    stats_exclude = NULL,
    marker_levels = NULL,
    legend_position = c("bottom","right","left","top","none"),
    save_plots = FALSE,
    save_dir = "",
    ncol = 3,
    sample_labeling = c("map_id", "sample_id_rep", "sample_id", "replicate"),
    engine = c("ggplot", "plotly")
){

  data <- .resolve_stats_summary_input(
    data = data,
    epk = epk,
    stats_summary = stats_summary,
    fn_name = "plot_qc_stats"
  )

  # Setup arguments
  engine <- match.arg(engine)
  legend_position <- match.arg(legend_position)
  sample_labeling <- match.arg(sample_labeling)
  show_legend <- legend_position != "none"
  should_save_plots <- isTRUE(save_plots) && engine == "ggplot"


  if (save_plots && engine == "plotly") {
    warning("save_plots is ignored for engine = 'plotly'. Render the widget in HTML or save using htmlwidgets::saveWidget().")
  }

  out <- list()

  if (should_save_plots && save_dir == "") {
    stop("Please provide a directory to save the plots.")
  }

  if (is.null(marker_levels)) marker_levels <- unique(data$marker)
  data$marker <- factor(data$marker, levels = marker_levels)

  if (!sample_labeling %in% names(data)) {
    stop("Column '", sample_labeling, "' was not found in `data`.")
  }

  if (is.null(stats)) {
    stats <- names(dplyr::select(data, where(is.numeric)))
  } else {
    stats <- intersect(stats, names(data))
    if (length(stats) == 0) stop("None of the requested stats exist in `data`.")
  }

  if (!is.null(stats_exclude)) {
    stats_exclude <- intersect(as.character(stats_exclude), names(data))
    stats <- setdiff(stats, stats_exclude)
    if (length(stats) == 0) {
      stop("No stats remain after applying `stats_exclude`.")
    }
  }

  if (!"condition" %in% names(data)) data$condition <- "all_conditions"

  if (is.null(condition)) {
    conditions <- unique(data$condition)
  } else if (identical(condition, "All")) {
    data$condition <- "All"
    conditions <- "All"
  } else {
    conditions <- condition
  }

  for (cond in conditions) {
    dfc <- dplyr::filter(data, .data$condition == cond)

    plot_list <- list()

    for (sts in stats) {
      p <- ggplot2::ggplot(dfc, ggplot2::aes(x = marker, y = .data[[sts]])) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = marker), outlier.shape = NA) +
        ggplot2::geom_point(
          ggplot2::aes(color = .data[[sample_labeling]]),
          position = ggplot2::position_jitter(width = 0.2),
          size = 2, alpha = 0.7
        ) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = gsub("_", " ", sts), x = NULL, y = NULL) +
        ggplot2::theme(
          plot.title = ggplot2::element_blank(),
          # IMPORTANT CHANGE:
          # keep legend ON for plotly, allow collecting for ggplot/patchwork
          legend.position = if (engine == "plotly" && show_legend) "right" else "none",
          axis.text.x  = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        ) +
        ggplot2::guides(fill = "none")

      plot_list[[sts]] <- p
    }

    # bottom-row x-axis labels (works for both engines)
    n <- length(plot_list)
    bottom_idx <- seq.int(from = max(1, n - ncol + 1), to = n)
    plot_list[bottom_idx] <- lapply(
      plot_list[bottom_idx],
      \(p) p + ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = ggplot2::element_line()
      )
    )
    plot_titles <- names(plot_list)

    if (engine == "ggplot") {
      composed <- patchwork::wrap_plots(plot_list, ncol = ncol) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = if (show_legend) legend_position else "none")

      out[[cond]] <- composed

      if (should_save_plots) {
        ggplot2::ggsave(
          filename = file.path(save_dir, paste0("QC_", cond, ".png")),
          plot = composed, width = 12, height = 7
        )
      }

    } else {
      # plotly engine: convert each ggplot and show legend only once
      p_list <- lapply(plot_list, function(p){
        pp <- plotly::ggplotly(p)
        # remove ggplot title and annotations (we'll add custom ones later)
        pp$x$layout$title <- NULL
        pp$x$layout$annotations <- NULL
        pp
      })

      p_list <- Map(function(p, i) {
        plotly::layout(p, showlegend = (show_legend && i == 1))
      }, p_list, seq_along(p_list))

      # compose subplots
      composed <- plotly::subplot(
        p_list,
        nrows = ceiling(length(p_list) / ncol),
        shareX = FALSE,
        shareY = FALSE
      )

      # turn off outlier/point display for plotly box traces
      composed <- plotly::style(
        composed,
        boxpoints = FALSE,
        traces = which(vapply(composed$x$data, function(tr) identical(tr$type, "box"), logical(1)))
      )

      # Give top margin space for titles
      composed <- plotly::layout(composed, margin = list(t = 60))

      n <- length(p_list)

      # helper to get layout key (xaxis, xaxis2, ...)
      ax_key <- function(prefix, i) if (i == 1) prefix else paste0(prefix, i)

      annotations <- lapply(seq_len(n), function(i) {
        xa <- composed$x$layout[[ax_key("xaxis", i)]]
        ya <- composed$x$layout[[ax_key("yaxis", i)]]

        xmid <- mean(unlist(xa$domain))
        ytop <- unlist(ya$domain)[2]

        list(
          text = plot_titles[[i]],
          x = xmid, #((i - 1) %% ncol) / ncol + 0.5 / ncol,
          y =ytop + 0.02, # 1 - (floor((i - 1) / ncol) / nrows) - 0.05,
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "bottom",
          showarrow = FALSE,
          font = list(size = 12)
        )
      })
      composed <- plotly::layout(
        composed,
        annotations = annotations,
        title = NULL
      )

      if(show_legend){
        composed <- plotly::layout(
          composed,
          legend = list(
            orientation = if (legend_position == "bottom") "h" else "v",
            x = 0,
            y = if (legend_position == "bottom") -0.08 else 1
          )
        )

      }else{
        composed <- plotly::layout(
          composed,
          showlegend = FALSE
        )
      }
      out[[cond]] <- composed
    }
  }

  if (length(out) == 1) return(out[[1]])
  return(out)
}


#' Plot scaling factors (msr) across samples
#'
#' Creates bar plots of scaling factors (\code{msr}) across samples, faceted by
#' marker. This function follows the same input framework as
#' \code{plot_qc_stats()} and supports \code{epk}, \code{stats_summary}, or
#' legacy \code{data} input.
#'
#' @param data Backward-compatible alias for \code{stats_summary}. A data frame
#'   containing at least \code{marker}, \code{sample_id_rep}, and \code{msr}.
#' @param epk Optional EPK object. If supplied,\code{epk$tables$stats_summary}
#'   is used as plotting input.
#' @param stats_summary Optional data frame containing scaling information. If
#'   provided, it takes precedence over \code{epk} and \code{data}.
#' @param condition Character vector specifying which conditions to plot.
#'   \itemize{
#'     \item \code{NULL}: plot all available conditions as-is (default)
#'     \item \code{"All"}: collapse all rows into a single condition named \code{"All"}
#'     \item otherwise: one or more condition names to subset
#'   }
#' @param marker_levels Optional character vector defining marker order for
#'   facet display/fill legend. If \code{NULL}, uses order in data.
#' @param markers_to_exclude Character vector of markers to remove before
#'   plotting. Default: \code{"INPUT"}.
#' @param legend_position Character string passed to
#'   \code{theme(legend.position = ...)}.
#' @param save_plots Logical; if \code{TRUE}, save PNG files.
#' @param save_dir Character path to directory where plots are saved when
#'   \code{save_plots = TRUE}.
#' @param ncol Integer; number of facet columns.
#' @param sample_labeling Character; column used for x-axis labels. One of
#'   \code{"sample_id_rep"}, \code{"sample_id"}, or \code{"map_id"}.
#'
#' @return If one condition is selected, returns one ggplot object. If multiple
#'   conditions are selected, returns a named list of ggplot objects.
#'
#' @examples
#' \dontrun{
#' p <- scaling_plot(
#'   epk = epk,
#'   condition = "All",
#'   legend_position = "none",
#'   ncol = 3,
#'   sample_labeling = "sample_id_rep"
#' )
#' p
#' }
#'
#' @importFrom ggplot2 ggplot geom_bar theme_bw labs geom_hline facet_wrap theme
#' @importFrom ggplot2 element_text ggsave aes
#' @importFrom dplyr filter
#' @export
scaling_plot <- function(
  data = NULL,
  epk = NULL,
  stats_summary = NULL,
  condition = NULL,
  marker_levels = NULL,
  markers_to_exclude = c("INPUT"),
  legend_position = c("bottom", "right", "left", "top", "none"),
  save_plots = FALSE,
  save_dir = "",
  ncol = 3,
  sample_labeling = c("sample_id_rep", "sample_id", "map_id")
) {

  data <- .resolve_stats_summary_input(
    data = data,
    epk = epk,
    stats_summary = stats_summary,
    fn_name = "scaling_plot"
  )

  legend_position <- match.arg(legend_position)
  sample_labeling <- match.arg(sample_labeling)

  required_cols <- c("marker", "msr", sample_labeling)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  if (!"condition" %in% names(data)) {
    data$condition <- "all_conditions"
  }

  if (!is.null(markers_to_exclude) && length(markers_to_exclude) > 0) {
    data <- dplyr::filter(data, !(.data$marker %in% markers_to_exclude))
  }

  if (nrow(data) == 0) {
    stop("No rows remain after filtering markers. Check 'markers_to_exclude'.")
  }

  if (is.null(marker_levels)) {
    marker_levels <- unique(data$marker)
  }
  data$marker <- factor(data$marker, levels = marker_levels)

  if (is.null(condition)) {
    conditions <- unique(data$condition)
  } else if (identical(condition, "All")) {
    data$condition <- "All"
    conditions <- "All"
  } else {
    conditions <- condition
  }

  out <- list()

  if (isTRUE(save_plots) && identical(save_dir, "")) {
    stop("Please provide a directory to save the plots.")
  }

  for (cond in conditions) {
    dfc <- dplyr::filter(data, .data$condition == cond)

    if (nrow(dfc) == 0) {
      warning("Condition '", cond, "' has no rows; skipping.")
      next
    }

    p <- ggplot2::ggplot(
      dfc,
      ggplot2::aes(x = .data[[sample_labeling]], y = .data$msr, fill = marker)
    ) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Scaling factors (msr) across samples",
        x = sample_labeling,
        y = "msr"
      ) +
      ggplot2::geom_hline(yintercept = 1) +
      ggplot2::facet_wrap(~marker, scales = "free_y", ncol = ncol) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = legend_position
      )

    out[[cond]] <- p

    if (isTRUE(save_plots)) {
      ggplot2::ggsave(
        filename = file.path(save_dir, paste0("scaling_", cond, ".png")),
        plot = p,
        width = 12,
        height = 7
      )
    }
  }

  if (length(out) == 0) {
    stop("No plots were created for the requested condition(s).")
  }
  if (length(out) == 1) return(out[[1]])
  out
}


#' Download ChromHMM coreMarks annotation files
#'
#' Downloads ChromHMM segmentation BED files from the Roadmap Epigenomics
#' repository (coreMarks joint model) if they are not already present locally.
#' Compressed files are automatically decompressed.
#'
#' @param annotations Character vector of ChromHMM BED filenames
#'   (e.g. "E107_15_coreMarks_hg38lift_mnemonics.bed").
#' @param dest_dir Destination directory where files will be stored.
#'   Defaults to "~/Epigenica/Data/chromHmm_annotation_files/".
#' @param base_url Base URL of the ChromHMM annotations repository.
#' @param overwrite Logical; whether to overwrite existing files.
#'
#' @return Invisibly returns a character vector of downloaded files.
#'
#' @examples
#' \dontrun{
#' download_chromhmm_annotations(
#'   annotations = c("E107_15_coreMarks_hg38lift_mnemonics.bed")
#' )
#' }
#'
#' @export
download_chromhmm_annotations <- function(
    annotations,
    dest_dir = "~/Epigenica/Data/chromHmm_annotation_files/",
    base_url = "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/",
    overwrite = FALSE
) {

  dest_dir <- path.expand(dest_dir)
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  local_files <- list.files(dest_dir, pattern = "\\.bed$", full.names = FALSE)
  missing_files <- setdiff(annotations, local_files)

  if (length(missing_files) == 0) {
    message("All ChromHMM annotations already present.")
    return(invisible(character()))
  }

  downloaded <- character()

  for (chm in missing_files) {

    gz_file  <- paste0(chm, ".gz")
    gz_path  <- file.path(dest_dir, gz_file)
    bed_path <- file.path(dest_dir, chm)

    message("Downloading ChromHMM annotation: ", chm)

    download.file(
      url      = paste0(base_url, gz_file),
      destfile = gz_path,
      mode     = "wb"
    )

    R.utils::gunzip(
      gz_path,
      destname = bed_path,
      overwrite = overwrite,
      remove = TRUE
    )

    downloaded <- c(downloaded, bed_path)
  }

  invisible(downloaded)
}



#' Create metadata dataframe from BigWig filenames
#' @param map_id_vector A character vector of mapping IDs (optional).
#' @param bw_files A character vector of BigWig file paths (optional).
#' @return A tibble containing metadata extracted from the filenames.
#' @examples
#' \dontrun{
#' bw_files <- c(
#'  "/path/to/project1_batch1_H3K4me3_rerun_sample1_rep1.hg38.scaled.bw",
#'  "/path/to/project1_batch1_H3K9me3_rerun_sample2_pooled.hg38.unscaled.bw"
#'  )
#'  create_metadata_df(bw_files = bw_files)
#'  }
#'  #  # A tibble: 2 × 8
#'  #  bw_file                                      project_id batch   marker   rerun_id sample_id
#'  #  <chr>                                        <chr>      <chr>   <chr>    <chr>     <chr>
#'  #1 /path/to/project1_batch1_H3K4me3_rerun_sa… project1    batch1  H3K4me3  rerun     sample1
#'  #2 /path/to/project1_batch1_H3K9me3_rerun_sa… project1    batch1  H3K9me3  rerun     sample2
#'  #  replicate genome  scaling
#'  #  <chr>     <chr>   <chr>
#'  #1 rep1      hg38    scaled
#'  #2 pooled    hg38    unscaled
#'
#' @importFrom stringr str_match
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @export
create_metadata_df <- function(map_id_vector= NULL, bw_files = NULL) {
  if(is.null(map_id_vector)){
    if(is.null(bw_files)){
      stop("Please provide either 'map_id_vector' or 'bw_files'.")
    } else {
      map_id_vector <- basename(bw_files)
      m <- stringr::str_match(
        map_id_vector,
        "^([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_(rep[0-9]+|pooled)\\.(hg38)\\.(scaled|unscaled)\\.bw$")
      # m columns: [1]=fullmatch, then capture groups 1..8
      bw_df <- tibble::tibble(
        bw_file    = map_id_vector,
        project_id = m[,2],
        batch      = m[,3],
        marker     = m[,4],
        rerun_id   = m[,5],
        sample_id  = m[,6],
        replicate  = m[,7],
        genome     = m[,8],
        scaling    = m[,9]
      )

    }
  }else{
    m <- stringr::str_match(
      map_id_vector,
      "^([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_(rep[0-9]+|pooled)\\.(hg38)$")
    # m columns: [1]=fullmatch, then capture groups 1..7
    bw_df <- tibble::tibble(
      map_id     = map_id_vector,
      project_id = m[,2],
      batch      = m[,3],
      marker     = m[,4],
      rerun_id   = m[,5],
      sample_id  = m[,6],
      replicate  = m[,7],
      genome     = m[,8]
    )
  }



  # Optional: flag non-matching files
  bw_df <- dplyr::mutate(
    bw_df,
    matched = !is.na(project_id)
  )

  return(bw_df)
}


#' Ensure a GTF annotation and derived BED files exist
#'
#' Downloads (and gunzips) a GTF file if it is missing (or if
#' \code{force_gtf_redownload = TRUE}), then creates one or both derived BED files
#' if they do not exist:
#' \itemize{
#'   \item \strong{Gene BED}: \code{chrom, start(0-based), end, gene_name, gene_id, score, strand, biotype}
#'   \item \strong{TSS ±2kb BED}: strand-aware TSS window per gene (\code{gene_name, gene_id} kept)
#' }
#'
#' The GTF is parsed using only rows with feature type \code{"gene"} (column V3).
#' TSS is defined as:
#' \itemize{
#'   \item \code{+} strand: \code{start} (0-based BED start)
#'   \item \code{-} strand: \code{end} (BED end)
#' }
#'
#' @param gtf_file Character. Path to the uncompressed GTF file.
#' @param gtf_url Character. URL to a gzipped GTF (\code{.gtf.gz}) to download when needed.
#' @param genes_bed Character. Path to the output gene BED file.
#' @param tss2k_bed Character. Path to the output TSS ±2kb BED file.
#' @param force_gtf_redownload Logical. If \code{TRUE}, re-downloads and overwrites \code{gtf_file}
#'   even if it already exists. Default: \code{FALSE}.
#' @param overwrite_beds Logical. If \code{TRUE}, (re)writes BED outputs even if they exist.
#'   Default: \code{FALSE}.
#'
#' @return An (invisible) named list with paths:
#' \code{list(gtf = gtf_file, genes_bed = genes_bed, tss2k_bed = tss2k_bed)}.
#'
#' @details
#' This function intentionally does \emph{not} re-download the GTF just because a BED
#' is missing (unless \code{force_gtf_redownload = TRUE}). If a BED is missing, it is
#' regenerated from the local GTF to ensure reproducibility and consistency.
#'
#' The \code{biotype} field is extracted from the attributes column using the key
#' \code{gene_type}. If your GTF uses a different key (e.g. \code{gene_biotype}),
#' adjust the regex accordingly.
#'
#' @examples
#' \dontrun{
#' paths <- ensure_gtf_and_beds(
#'   gtf_file  = "~/Epigenica/Data/gencode.v38.annotation.gtf",
#'   gtf_url   = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz",
#'   genes_bed = "~/Epigenica/Data/genes.hg38.bed",
#'   tss2k_bed = "~/Epigenica/Data/genes_tss_2kb.hg38.bed"
#' )
#' genes_coord   <- paths$genes_bed
#' tss_2kb_coord <- paths$tss2k_bed
#' }
#'
#' @importFrom utils download.file
#' @importFrom stringr str_match
#' @export
ensure_gtf_and_beds <- function(
    gtf_file,
    gtf_url,
    genes_bed,
    tss2k_bed,
    force_gtf_redownload = FALSE,
    overwrite_beds = FALSE
) {
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required. Install it with install.packages('stringr').")
  }
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package 'R.utils' is required. Install it with install.packages('R.utils').")
  }

  dir.create(dirname(gtf_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(genes_bed), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(tss2k_bed), recursive = TRUE, showWarnings = FALSE)

  # ---- 1) Ensure GTF ----
  if (force_gtf_redownload || !file.exists(gtf_file)) {
    gz <- paste0(gtf_file, ".gz")
    utils::download.file(url = gtf_url, destfile = gz, mode = "wb")
    R.utils::gunzip(gz, destname = gtf_file, overwrite = TRUE, remove = TRUE)
  }

  # ---- 2) Decide what to build ----
  need_genes <- overwrite_beds || !file.exists(genes_bed)
  need_tss2k <- overwrite_beds || !file.exists(tss2k_bed)

  if (!need_genes && !need_tss2k) {
    return(invisible(list(gtf = gtf_file, genes_bed = genes_bed, tss2k_bed = tss2k_bed)))
  }

  # ---- 3) Read and parse once ----
  gtf_df <- read.table(
    gtf_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "#",
    quote = ""
  )

  # Keep only 'gene' features
  gtf_df <- gtf_df[gtf_df$V3 == "gene", , drop = FALSE]

  gene_name <- stringr::str_match(gtf_df$V9, "gene_name\\s+([^;]+);")[, 2]
  gene_id   <- stringr::str_match(gtf_df$V9, "gene_id\\s+([^;]+);")[, 2]
  biotype   <- stringr::str_match(gtf_df$V9, "gene_type\\s+([^;]+);")[, 2]

  bed_df <- data.frame(
    chrom     = gtf_df$V1,
    start     = gtf_df$V4 - 1,  # convert to 0-based BED start
    end       = gtf_df$V5,
    gene_name = gene_name,
    gene_id   = gene_id,
    score     = 0,
    strand    = gtf_df$V7,
    biotype   = biotype,
    stringsAsFactors = FALSE
  )

  # ---- 4) Write outputs (as needed) ----
  if (need_genes) {
    write.table(
      bed_df,
      file = genes_bed,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }

  if (need_tss2k) {
    # Strand-aware TSS: + => start; - => end
    tss <- ifelse(bed_df$strand == "+", bed_df$start, bed_df$end)

    tss_bed <- data.frame(
      chrom     = bed_df$chrom,
      start     = pmax(0, tss - 2000),
      end       = tss + 2000,
      gene_name = bed_df$gene_name,
      gene_id   = bed_df$gene_id,
      score     = 0,
      strand    = bed_df$strand,
      stringsAsFactors = FALSE
    )

    write.table(
      tss_bed,
      file = tss2k_bed,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }

  invisible(list(gtf = gtf_file, genes_bed = genes_bed, tss2k_bed = tss2k_bed))
}

#' Interactive enrichment profile plot
#'
#' Generates an interactive multi-sample enrichment profile plot from binned
#' signal values. Lines are colored by sample class (control vs non-control),
#' and hovered samples are highlighted while non-hovered samples are dimmed.
#'
#' @param df A data frame containing at least the columns 
#'   \\code{index}, \\code{mean}, and \\code{sample}.
#' @param bin_size Numeric; bin size in base pairs. Default: \\code{100}.
#' @param window_bp Numeric; half window size in base pairs. Default:
#'   \\code{2500} (i.e. +/- 2.5 kb).
#' @param highlight_opacity Numeric in [0,1]; opacity applied to non-hovered
#'   traces. Default: \\code{0.5}.
#' @param mid_coord Character; reference point used for center label and x-axis
#'   title. One of \\code{"center"} or \\code{"start"}. Default:
#'   \\code{"center"}.
#' @param plot_title Character; title for the plot. Default:
#'   \\code{"Enrichment Profile"}.
#' @param control_samples Optional character vector of sample names to treat as
#'   controls. If provided, takes precedence over \\code{control_pattern}.
#' @param control_pattern Optional regex pattern used to classify controls when
#'   \\code{control_samples} is \\code{NULL}.
#' @param control_color Character color for control traces. Default:
#'   \\code{"rgba(0,0,0,0.90)"}.
#' @param other_color Character color for non-control traces. Default:
#'   \\code{"rgba(120,180,255,0.55)"}.
#' @param hover_color Character color used for highlighted trace on hover.
#'   Default: \\code{"rgba(0,120,255,1)"}.
#' @param group_by Optional character scalar naming a column used for color
#'   grouping (e.g. \\code{"condition"}, \\code{"batch"}). If provided,
#'   samples are colored by this grouping instead of control vs non-control.
#' @param metadata Optional data frame containing sample annotations to join
#'   into \\code{df} before plotting (for example, condition assignments).
#' @param metadata_sample_col Character scalar naming the sample-id column in
#'   \\code{metadata}. Default: \\code{"sample"}.
#' @param group_palette Optional named character vector mapping group names to
#'   colors. When named, names should match the observed \\code{group_by}
#'   levels. Missing levels are auto-filled (with warning); extra names are
#'   ignored (with warning).
#' @param control_overrides_group Logical; when \\code{TRUE} and controls are
#'   defined, control samples use \\code{control_color} even when
#'   \\code{group_by} is provided. Default: \\code{FALSE}.
#'
#' @return A \\code{plotly} htmlwidget.
#'
#' @details
#' The x-axis is computed as genomic distance in bp from the mid-window and
#' labeled as \\code{-window}, \\code{mid_coord}, and \\code{+window}. A dashed
#' vertical line marks \\code{0 bp}. Hover highlighting remains sample-specific,
#' even when base colors are assigned by group.
#'
#' @examples
#' \\dontrun{
#' set.seed(1)
#' ex <- expand.grid(index = 1:51, sample = c("INPUT_A1", "H3K4me3_A1"))
#' ex$mean <- c(dnorm(seq(-2.5, 2.5, length.out = 51), sd = 0.8),
#'              dnorm(seq(-2.5, 2.5, length.out = 51), sd = 0.5))
#'
#' p <- plot_enrichment_interactive(
#'   ex,
#'   bin_size = 100,
#'   window_bp = 2500,
#'   mid_coord = "start",
#'   group_by = "condition",
#'   group_palette = c(INPUT = "black", H3K4me3 = "#1f77b4"),
#'   control_pattern = "^INPUT",
#'   control_overrides_group = TRUE
#' )
#' p
#' }
#'
#' @export
plot_enrichment_interactive <- function(
    df,
    bin_size = 100,
    window_bp = 2500,
    highlight_opacity = 0.5,
    mid_coord = c("center", "start"),
    plot_title = "Enrichment Profile",
    control_samples = NULL,
    control_pattern = NULL,
    control_color = "rgba(0,0,0,0.90)",
    other_color = "rgba(120,180,255,0.55)",
    hover_color = "rgba(0,120,255,1)",
    group_by = NULL,
    metadata = NULL,
    metadata_sample_col = "sample",
    group_palette = NULL,
    control_overrides_group = FALSE
) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame. Got: ", paste(class(df), collapse = "/"))
  }

  required_cols <- c("index", "mean", "sample")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("'df' is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(bin_size) || length(bin_size) != 1 || is.na(bin_size) || bin_size <= 0) {
    stop("'bin_size' must be a single positive numeric value.")
  }
  if (!is.numeric(window_bp) || length(window_bp) != 1 || is.na(window_bp) || window_bp <= 0) {
    stop("'window_bp' must be a single positive numeric value.")
  }
  if (!is.numeric(highlight_opacity) || length(highlight_opacity) != 1 ||
      is.na(highlight_opacity) || highlight_opacity < 0 || highlight_opacity > 1) {
    stop("'highlight_opacity' must be a single numeric value in [0, 1].")
  }
  if (!is.null(group_by) && (!is.character(group_by) || length(group_by) != 1)) {
    stop("'group_by' must be NULL or a single column name.")
  }
  if (!is.character(metadata_sample_col) || length(metadata_sample_col) != 1) {
    stop("'metadata_sample_col' must be a single column name.")
  }
  if (!is.logical(control_overrides_group) || length(control_overrides_group) != 1) {
    stop("'control_overrides_group' must be TRUE or FALSE.")
  }

  mid_coord <- match.arg(mid_coord)

  .normalize_plotly_color <- function(x) {
    x <- as.character(x)
    out <- x
    is_rgba <- grepl("^rgba\\(", x)

    if (any(is_rgba)) {
      parsed <- regexec(
        "^rgba\\((\\d+),(\\d+),(\\d+),([0-9]*\\.?[0-9]+)\\)$",
        gsub("\\s+", "", x[is_rgba])
      )
      m <- regmatches(gsub("\\s+", "", x[is_rgba]), parsed)

      out_rgba <- vapply(m, function(v) {
        if (length(v) != 5) return(NA_character_)
        r <- as.numeric(v[2])
        g <- as.numeric(v[3])
        b <- as.numeric(v[4])
        a <- as.numeric(v[5])
        if (any(is.na(c(r, g, b, a)))) return(NA_character_)
        grDevices::rgb(r, g, b, alpha = a * 255, maxColorValue = 255)
      }, character(1))

      out[is_rgba] <- out_rgba
    }

    out
  }

  control_color <- .normalize_plotly_color(control_color)
  other_color <- .normalize_plotly_color(other_color)
  hover_color <- .normalize_plotly_color(hover_color)

  df2 <- dplyr::mutate(
    df,
    index = suppressWarnings(as.numeric(index)),
    mean = suppressWarnings(as.numeric(mean)),
    sample = trimws(as.character(sample))
  )

  if (anyNA(df2$index)) {
    stop("'index' contains non-numeric values (or NA after coercion).")
  }
  if (anyNA(df2$mean)) {
    stop("'mean' contains non-numeric values (or NA after coercion).")
  }
  if (any(df2$sample == "")) {
    stop("'sample' contains empty values.")
  }

  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) {
      stop("'metadata' must be a data.frame when provided.")
    }
    if (!metadata_sample_col %in% names(metadata)) {
      stop("'metadata' is missing join column: ", metadata_sample_col)
    }

    md <- as.data.frame(metadata, stringsAsFactors = FALSE)
    md[[metadata_sample_col]] <- trimws(as.character(md[[metadata_sample_col]]))

    if (any(md[[metadata_sample_col]] == "")) {
      stop("'metadata' contains empty sample identifiers in column '", metadata_sample_col, "'.")
    }
    if (anyDuplicated(md[[metadata_sample_col]]) > 0) {
      stop("'metadata' contains duplicated sample identifiers in column '", metadata_sample_col, "'.")
    }

    idx <- match(df2$sample, md[[metadata_sample_col]])
    md_cols <- setdiff(names(md), metadata_sample_col)
    for (cn in md_cols) {
      df2[[cn]] <- md[[cn]][idx]
    }
  }

  is_control_vec <- rep(FALSE, nrow(df2))
  if (!is.null(control_samples)) {
    is_control_vec <- df2$sample %in% as.character(control_samples)
  } else if (!is.null(control_pattern)) {
    is_control_vec <- grepl(control_pattern, df2$sample)
  }

  df2 <- dplyr::mutate(
    df2,
    is_control = is_control_vec,
    pos_bp = -window_bp + (index - 1) * bin_size,
    sample = factor(sample)
  )

  sample_info <- dplyr::distinct(df2, sample, is_control)

  if (!is.null(group_by)) {
    if (!group_by %in% names(df2)) {
      stop("'group_by' column not found in data: ", group_by)
    }

    sample_group <- dplyr::distinct(df2, sample, group = .data[[group_by]])
    sample_group$group <- as.character(sample_group$group)
    sample_group$group[is.na(sample_group$group) | trimws(sample_group$group) == ""] <- "Unknown"

    if (any(duplicated(sample_group$sample))) {
      stop("Each sample must map to exactly one group in column '", group_by, "'.")
    }

    grp_levels <- sort(unique(sample_group$group))

    if (is.null(group_palette)) {
      grp_colors <- grDevices::hcl.colors(length(grp_levels), palette = "Set 2")
      group_color_map <- stats::setNames(grp_colors, grp_levels)
    } else {
      if (!is.character(group_palette)) {
        stop("'group_palette' must be a character vector of colors.")
      }
      group_palette <- .normalize_plotly_color(group_palette)
      if (is.null(names(group_palette))) {
        if (length(group_palette) < length(grp_levels)) {
          stop("Unnamed 'group_palette' must have at least as many colors as unique groups.")
        }
        group_color_map <- stats::setNames(group_palette[seq_along(grp_levels)], grp_levels)
      } else {
        group_color_map <- group_palette
        missing_groups <- setdiff(grp_levels, names(group_color_map))
        extra_groups <- setdiff(names(group_color_map), grp_levels)

        if (length(extra_groups) > 0) {
          warning(
            "'group_palette' contains names not present in observed groups: ",
            paste(extra_groups, collapse = ", "),
            ". They will be ignored."
          )
          group_color_map <- group_color_map[setdiff(names(group_color_map), extra_groups)]
        }

        if (length(missing_groups) > 0) {
          warning(
            "'group_palette' is missing colors for group level(s): ",
            paste(missing_groups, collapse = ", "),
            ". Auto-generating colors for missing levels."
          )
          auto_cols <- grDevices::hcl.colors(length(missing_groups), palette = "Set 2")
          group_color_map <- c(group_color_map, stats::setNames(auto_cols, missing_groups))
        }
      }
    }

    sample_color_map <- stats::setNames(
      group_color_map[sample_group$group],
      as.character(sample_group$sample)
    )

    sample_info <- dplyr::left_join(sample_info, sample_group, by = "sample")

    if (isTRUE(control_overrides_group)) {
      sample_color_map[as.character(sample_info$sample[sample_info$is_control])] <- control_color
    }

    color_map <- sample_color_map
  } else {
    color_map <- stats::setNames(
      ifelse(sample_info$is_control, control_color, other_color),
      as.character(sample_info$sample)
    )
  }

  dfk <- plotly::highlight_key(df2, ~sample)

  window_kb <- formatC(window_bp / 1000, format = "fg", digits = 3)
  x_tick_vals <- c(-window_bp, 0, window_bp)
  x_tick_text <- c(
    paste0("-", window_kb, " kb"),
    mid_coord,
    paste0("+", window_kb, " kb")
  )

  hover_template <- "<b>%{fullData.name}</b><br>mean: %{y:.3f}<extra></extra>"
  if (!is.null(group_by)) {
    hover_template <- paste0(
      "<b>%{fullData.name}</b><br>",
      group_by,
      ": %{customdata}<br>mean: %{y:.3f}<extra></extra>"
    )
  }

  if (!is.null(group_by)) {
    df2$..group_plotly <- as.character(df2[[group_by]])
    df2$..group_plotly[is.na(df2$..group_plotly) | trimws(df2$..group_plotly) == ""] <- "Unknown"
    dfk <- plotly::highlight_key(df2, ~sample)
  }

  if (!is.null(group_by)) {
    p <- plotly::plot_ly(
      data = dfk,
      x = ~pos_bp,
      y = ~mean,
      type = "scatter",
      mode = "lines",
      color = ~sample,
      colors = color_map,
      line = list(width = 1),
      customdata = ~..group_plotly,
      hovertemplate = hover_template
    )
  } else {
    p <- plotly::plot_ly(
      data = dfk,
      x = ~pos_bp,
      y = ~mean,
      type = "scatter",
      mode = "lines",
      color = ~sample,
      colors = color_map,
      line = list(width = 1),
      hovertemplate = hover_template
    )
  }

  p <- p %>%
    plotly::layout(
      showlegend = FALSE,
      title = plot_title,
      xaxis = list(
        title = paste0("Distance from ", mid_coord, " (bp)"),
        tickmode = "array",
        tickvals = x_tick_vals,
        ticktext = x_tick_text
      ),
      yaxis = list(title = "Mean signal"),
      shapes = list(list(
        type = "line",
        x0 = 0,
        x1 = 0,
        y0 = 0,
        y1 = 1,
        yref = "paper",
        line = list(color = "rgba(100,100,100,0.8)", dash = "dash", width = 1)
      ))
    ) %>%
    plotly::highlight(
      on = "plotly_hover",
      off = "plotly_doubleclick",
      persistent = FALSE,
      opacityDim = highlight_opacity,
      selected = list(line = list(color = hover_color, width = 3))
    )

  p
}

#' Print an EPK object
#'
#' Prints a compact summary of an \code{EPK} object without printing large tables.
#'
#' @param x An \code{EPK} object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
#' @method print EPK
print.EPK <- function(x, ...) {
  cat("EPK object\n")
  cat("-----------\n")

  # Extract and display markers
  markers <- NULL
  if (!is.null(x$tables$stats_summary) && "marker" %in% names(x$tables$stats_summary)) {
    markers <- unique(x$tables$stats_summary$marker)
    markers <- markers[!is.na(markers) & markers != ""]
  } else if (!is.null(x$mse)) {
    exps <- MultiAssayExperiment::experiments(x$mse)
    if (length(exps) > 0) {
      se <- exps[[1]]
      markers <- names(SummarizedExperiment::assays(se))
    }
  }
  if (!is.null(markers) && length(markers) > 0) {
    cat("* Markers:", paste(markers, collapse = ", "), "\n")
  }

  # MultiAssayExperiment summary
  if (!is.null(x$mse)) {
    exps <- MultiAssayExperiment::experiments(x$mse)
    cat("* Experiments:\n")
    for (nm in names(exps)) {
      se <- exps[[nm]]
      cat(sprintf(
        "  - %s: %d features × %d samples, %d assays\n",
        nm, nrow(se), ncol(se), length(SummarizedExperiment::assays(se))
      ))
    }
  }

  # Tables: only dimensions, not contents
  if (!is.null(x$tables) && length(x$tables) > 0) {
    cat("* Tables:\n")
    for (nm in names(x$tables)) {
      obj <- x$tables[[nm]]
      if (is.data.frame(obj)) {
        cat(sprintf("  - %s: %d rows × %d cols\n", nm, nrow(obj), ncol(obj)))
      } else if (!is.null(obj)) {
        cat(sprintf("  - %s: <%s>\n", nm, class(obj)[1]))
      }
    }
  }


  # Derived info (print before annotation to ensure visibility)
  if (!is.null(x$derived) && length(x$derived) > 0) {
    cat("* Derived:")
    for (nm in names(x$derived)) {
      obj <- x$derived[[nm]]
      if (is.list(obj)) {
        cat(sprintf("\n  - %s: list with %d entries", nm, length(obj)))
      } else if (!is.null(obj)) {
        cat(sprintf("\n  - %s: <%s>", nm, class(obj)[1]))
      }
    }
    cat("\n")
  }

  # Annotation: features annotation
  if (!is.null(x$annotation) && length(x$annotation) > 0) {
    cat("* Annotation:\n")
    for (nm in names(x$annotation)) {
      obj <- x$annotation[[nm]]
      if (is.data.frame(obj)) {
        cat(sprintf("  - %s: %d rows × %d cols\n", nm, nrow(obj), ncol(obj)))
      } else if (!is.null(obj)) {
        cat(sprintf("  - %s: <%s>\n", nm, class(obj)[1]))
      }
    }
  }

  # Other payloads
  if (!is.null(x$enrichment_results)) {
    cat(sprintf("* Enrichment results: %d tables\n", length(x$enrichment_results)))
  }
  if (!is.null(x$qc_plots)) {
    cat(sprintf("* QC plots: <%s>\n", class(x$qc_plots)[1]))
  }
  if (!is.null(x$enrichment_plots)) {
    cat(sprintf("* Enrichment plots: %d\n", length(x$enrichment_plots)))
  }
  if (!is.null(x$annotation)) {
    cat(sprintf("* Annotation GRanges: %d\n", length(x$annotation)))
  }

  if (!is.null(x$provenance$created)) {
    cat("* Created:", format(x$provenance$created), "\n")
  }

  invisible(x)
}

#' Compute sample–sample correlation for one marker
#'
#' Computes a correlation matrix between samples for a single assay
#' (i.e. one marker) stored in a \code{SummarizedExperiment}.
#'
#' @param se A \code{SummarizedExperiment}.
#' @param assay_name Name of the assay (marker).
#' @param method Correlation method. One of \code{"pearson"} or \code{"spearman"}.
#' @param transform Optional transformation applied before correlation.
#'   One of \code{"none"} or \code{"log1p"}.
#'
#' @return A numeric matrix (samples × samples) of correlations.
#' @export
compute_sample_cor <- function(epk, method = "pearson", transform = c("none", "log1p")) {
  transform <- match.arg(transform)
  exp_names <- names(MultiAssayExperiment::experiments(epk$mse))
  sample_cor_list <- setNames(
    lapply(exp_names, function(exp) {
      se <- MultiAssayExperiment::experiments(epk$mse)[[exp]]
      assay_names <- names(SummarizedExperiment::assays(se))
      setNames(
        lapply(assay_names, function(a) {
          m <- SummarizedExperiment::assay(se, a)
          if (transform == "log1p") m <- log1p(m)
          stats::cor(m, use = "pairwise.complete.obs", method = method)
        }),
        assay_names
      )
    }),
    exp_names
  )
  epk$derived <- c(epk$derived, list(sample_cor = sample_cor_list))
  epk
}

#' Compute sample–sample correlations for all assays in a specified experiment of an EPK object
#'
#' @param epk An EPK object (S3 list with an `mse` slot).
#' @param exp_name Name of the experiment (e.g. "protein_coding").
#' @param method Correlation method. One of "pearson" or "spearman".
#' @param transform Optional transformation applied before correlation.
#' @return The EPK object with a new slot epk$derived$all_cor containing the correlations for the specified experiment.
#' ## Example usage:
#' @examples
#' # Compute and add all correlations for a specific experiment
#' epk <- compute_all_cor(epk, exp_name = "genes", method = "pearson", transform = "none")
#' # Access correlations:
#' # epk$derived$all_cor[["genes"]]
#' @export
compute_all_cor <- function(epk, exp_name, method = "pearson", transform = c("none", "log1p")) {
  transform <- match.arg(transform)
  se <- MultiAssayExperiment::experiments(epk$mse)[[exp_name]]
  assay_names <- names(SummarizedExperiment::assays(se))
  all_cor_list <- setNames(
    lapply(assay_names, function(a) {
      m <- SummarizedExperiment::assay(se, a)
      if (transform == "log1p") m <- log1p(m)
      stats::cor(m, use = "pairwise.complete.obs", method = method)
    }),
    assay_names
  )
  if (is.null(epk$derived)) epk$derived <- list()
  if (is.null(epk$derived$all_cor)) epk$derived$all_cor <- list()
  epk$derived$all_cor[[exp_name]] <- all_cor_list
  epk
}



#' Calculate beta values from M2 and M3 assays
#'
#' Computes beta values as the ratio \eqn{(M2 + p) / (M2 + M3 + 2p)} where
#' \eqn{p} is the pseudocount.
#' Regions where both M2 and M3 are zero are set to \code{NA} to avoid bias
#' from the pseudocount.
#'
#' @param epk An EPK object (S3 list with an \code{mse} slot).
#' @param feature Character; name of the experiment in the
#'   \code{MultiAssayExperiment} (e.g. \code{"cpg_islands"}, \code{"union_peaks"}).
#' @param pseudocount Numeric; pseudocount added to numerator and denominator
#'   (default 1). The denominator receives \code{2 * pseudocount}.
#' @param method Character; one of \code{"raw"}, \code{"quantile"}, or \code{"both"}.
#'   \code{"quantile"} applies quantile normalisation to mbdseq and cxxc separately
#'   before computing the ratio (requires the \pkg{preprocessCore} package).
#'   \code{"both"} computes and stores both raw and quantile-normalized beta assays.
#' @param mbdseq Character; name of the MBD-seq assay (default \code{"mbdseq"}).
#' @param cxxc Character; name of the CXXC assay (default \code{"cxxc"}).
#'
#' @return The EPK object with a new assay added to the specified experiment:
#'   \code{"beta_raw"} or \code{"beta_qn"} depending on \code{method}.
#'
#' @examples
#' \dontrun{
#' epk <- readRDS("project.epk.rds")
#'
#' # Raw beta values on CpG islands (default)
#' epk <- calculate_beta_val(epk)
#'
#' # Quantile-normalised beta on union peaks
#' epk <- calculate_beta_val(epk, feature = "union_peaks", method = "quantile")
#'
#' # Custom pseudocount
#' epk <- calculate_beta_val(epk, pseudocount = 0.5, method = "raw")
#' }
#'
#' @export
calculate_beta_val <- function(epk,
                               feature = "cpg_islands",
                               pseudocount = 1,
                               method = c("raw", "quantile", "both"),
                               mbdseq = "mbdseq",
                               cxxc = "cxxc") {
  method <- match.arg(method)

  ## ---- validate EPK and extract SE ----
  if (!inherits(epk, "EPK") && is.null(epk$mse)) {
    stop("`epk` must be an EPK object with an `mse` slot.", call. = FALSE)
  }

  exp_names <- names(MultiAssayExperiment::experiments(epk$mse))
  if (!feature %in% exp_names) {
    stop(
      sprintf("Feature '%s' not found in EPK. Available: %s",
              feature, paste(exp_names, collapse = ", ")),
      call. = FALSE
    )
  }

  se <- MultiAssayExperiment::experiments(epk$mse)[[feature]]

  avail <- names(SummarizedExperiment::assays(se))
  missing <- setdiff(c(mbdseq, cxxc), avail)
  if (length(missing) > 0) {
    stop(
      sprintf("Assay(s) not found in '%s': %s. Available: %s",
              feature, paste(missing, collapse = ", "),
              paste(avail, collapse = ", ")),
      call. = FALSE
    )
  }

  ## ---- extract matrices ----
  m2 <- as.matrix(SummarizedExperiment::assay(se, mbdseq))
  m3 <- as.matrix(SummarizedExperiment::assay(se, cxxc))

  stopifnot(
    all(dim(m2) == dim(m3)),
    all(rownames(m2) == rownames(m3)),
    all(colnames(m2) == colnames(m3))
  )

  ## ---- mask double-zero regions ----
  both_zero <- m2 == 0 & m3 == 0
  m2[both_zero] <- NA
  m3[both_zero] <- NA

  if (method == "raw") {
      beta <- (m2 + pseudocount) / (m2 + m3 + 2 * pseudocount)
      assay_label <- paste0(mbdseq, "_", cxxc, "_raw")
      SummarizedExperiment::assay(se, assay_label) <- beta
      MultiAssayExperiment::experiments(epk$mse)[[feature]] <- se
      message(sprintf(
        "[beta] Added assay '%s' to experiment '%s' (%d features x %d samples).",
        assay_label, feature, nrow(beta), ncol(beta)
      ))
      return(epk)
  } else if (method == "quantile") {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop(
        "Package 'preprocessCore' is required for method = \"quantile\". ",
        "Install it with: BiocManager::install(\"preprocessCore\")",
        call. = FALSE
      )
    }
    m2_qn <- preprocessCore::normalize.quantiles(m2)
    m3_qn <- preprocessCore::normalize.quantiles(m3)
    dimnames(m2_qn) <- dimnames(m2)
    dimnames(m3_qn) <- dimnames(m3)
    beta <- (m2_qn + pseudocount) / (m2_qn + m3_qn + 2 * pseudocount)
    assay_label <- paste0(mbdseq, "_", cxxc, "_qn")
    SummarizedExperiment::assay(se, assay_label) <- beta
    MultiAssayExperiment::experiments(epk$mse)[[feature]] <- se
    message(sprintf(
      "[beta] Added assay '%s' to experiment '%s' (%d features x %d samples).",
      assay_label, feature, nrow(beta), ncol(beta)
    ))
    return(epk)
  } else if (method == "both") {
    # Raw
    beta_raw <- (m2 + pseudocount) / (m2 + m3 + 2 * pseudocount)
    assay_label_raw <- paste0(mbdseq, "_", cxxc, "_raw")
    SummarizedExperiment::assay(se, assay_label_raw) <- beta_raw
    # Quantile
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop(
        "Package 'preprocessCore' is required for method = \"both\" (quantile part). ",
        "Install it with: BiocManager::install(\"preprocessCore\")",
        call. = FALSE
      )
    }
    m2_qn <- preprocessCore::normalize.quantiles(m2)
    m3_qn <- preprocessCore::normalize.quantiles(m3)
    dimnames(m2_qn) <- dimnames(m2)
    dimnames(m3_qn) <- dimnames(m3)
    beta_qn <- (m2_qn + pseudocount) / (m2_qn + m3_qn + 2 * pseudocount)
    assay_label_qn <- paste0(mbdseq, "_", cxxc, "_qn")
    SummarizedExperiment::assay(se, assay_label_qn) <- beta_qn
    MultiAssayExperiment::experiments(epk$mse)[[feature]] <- se
    message(sprintf(
      "[beta] Added assays '%s' and '%s' to experiment '%s' (%d features x %d samples).",
      assay_label_raw, assay_label_qn, feature, nrow(beta_raw), ncol(beta_raw)
    ))
    return(epk)
  }

  ## ---- store back into EPK ----
  SummarizedExperiment::assay(se, assay_label) <- beta
  MultiAssayExperiment::experiments(epk$mse)[[feature]] <- se

  message(sprintf(
    "[beta] Added assay '%s' to experiment '%s' (%d features x %d samples).",
    assay_label, feature, nrow(beta), ncol(beta)
  ))

  epk
}

#' Plot a correlation heatmap for one or more markers in an EPK experiment
#'
#' @param epk An EPK object.
#' @param exp_name Name of the experiment (e.g. "protein_coding").
#' @param marker Name or vector of marker(s)/assay(s) (e.g. "H3K4me3" or c("5mC", "CXXC")).
#' @return A ComplexHeatmap::Heatmap object (single marker) or a list of Heatmap objects (multiple markers).
#' @examples
#' # Compute correlations first:
#' epk <- compute_all_cor(epk, exp_name = "genes")
#' # Plot heatmap for a single marker
#' heatmap_cor_marker(epk, exp_name = "genes", marker = "H3K4me3")
#' # Plot heatmaps for multiple markers
#' heatmap_cor_marker(epk, exp_name = "genes", marker = c("5mC", "CXXC"))
#' @export
heatmap_cor_marker <- function(epk, exp_name, marker) {
  # Helper for a single marker
  plot_one <- function(m) {
    cor_mat <- NULL
    if (!is.null(epk$derived) && !is.null(epk$derived$all_cor) &&
        !is.null(epk$derived$all_cor[[exp_name]]) &&
        !is.null(epk$derived$all_cor[[exp_name]][[m]])) {
      cor_mat <- epk$derived$all_cor[[exp_name]][[m]]
    } else {
      stop(sprintf(
        "No correlation matrix found for experiment '%s', marker '%s'.\nRun compute_all_cor(epk, '%s') first.",
        exp_name, m, exp_name))
    }
    stopifnot(is.matrix(cor_mat), !is.null(rownames(cor_mat)), !is.null(colnames(cor_mat)))
    ord <- hclust(as.dist(1 - cor_mat))$order
    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("white", "#fdae61", "#d73027"))
    ComplexHeatmap::Heatmap(
      cor_mat[ord, ord, drop = FALSE],
      name = paste0("cor_", m),
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      column_title = paste0("Sample–sample correlation, ", exp_name, ": ", m),
      heatmap_legend_param = list(title = "Correlation")
    )
  }
  if (length(marker) == 1) {
    plot_one(marker)
  } else {
    lapply(marker, plot_one)
  }
}

