#' Extract marker names from sample IDs
#' @param id A character string representing the sample ID.
#' @param markers A vector of marker names to search for within the sample ID.
#' @return A character vector containing the extracted marker names.
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
compute_sample_cor <- function(se, assay_name,
                               method = "pearson",
                               transform = c("none", "log1p")) {
  transform <- match.arg(transform)
  m <- SummarizedExperiment::assay(se, assay_name)
  if (transform == "log1p") m <- log1p(m)
  stats::cor(m, use = "pairwise.complete.obs", method = method)
}


#' Compute sample–sample correlations for all markers
#'
#' Computes sample–sample correlation matrices for all assays (markers)
#' within one experiment of a \code{MultiAssayExperiment}.
#'
#' @param mse A \code{MultiAssayExperiment}.
#' @param exp_name Name of the experiment (e.g. \code{"protein_coding"}).
#' @param method Correlation method. One of \code{"pearson"} or \code{"spearman"}.
#' @param transform Optional transformation applied before correlation.
#'
#' @return A named list of correlation matrices, one per marker.
#' @export
compute_all_cor <- function(mse, exp_name,
                            method = "pearson",
                            transform = "log1p") {
  se <- MultiAssayExperiment::experiments(mse)[[exp_name]]
  assay_names <- names(SummarizedExperiment::assays(se))
  setNames(
    lapply(assay_names, \(a) compute_sample_cor(se, a, method, transform)),
    assay_names
  )
}
