#' Create an EPK (Epigenomic Package) object from BigWig files and annotations
#'
#' High-level wrapper that simplifies EPK object creation. Accepts flexible input formats:
#' \itemize{
#'   \item \strong{Explicit mode}: \code{bw_files}, \code{annotations}, and \code{stats_summary}
#'   \item \strong{Path mode}: \code{pipeline_output_path} and \code{annotations} (auto-discovers files)
#' }
#'
#' @param bw_files Character vector of BigWig file paths. Required for explicit mode.
#'   Ignored if \code{pipeline_output_path} is provided.
#'
#' @param annotations Either a \code{GRanges} object, path to a BED file, or a data frame
#'   with columns \code{chr}, \code{start}, \code{end}, and optionally \code{feature_id}.
#'   Can also be a named list of such objects/paths for multiple annotation sets.
#'   Required in both modes.
#'
#' @param stats_summary A data frame containing QC statistics (e.g., from \code{toy_stats_summary}).
#'   Must include column \code{map_id} in format
#'   \code{(.*_(pooled|rep[0-9]+).genome)} so BigWig files can be checked against
#'   stats rows. Required for explicit mode.
#'   Ignored if \code{pipeline_output_path} is provided.
#'
#' @param pipeline_output_path Character path to pipeline output directory containing:
#'   \itemize{
#'     \item \code{minute_output/bigwig/} with \code{*.bw} files
#'     \item \code{minute_output/reports/stats_summary.txt}
#'   }
#'   If provided, \code{bw_files} and \code{stats_summary} are ignored.
#'
#' @param genome Character; genome version (default \code{"hg38"}).
#'   Used for validation and metadata.
#'
#' @param markers_to_exclude Character vector of marker names to skip
#'   (default \code{c("INPUT")}). Useful for excluding control tracks.
#'
#' @param experiment_names Character vector naming the experiments. If a single
#'   annotation is provided, defaults to \code{"protein_coding"}. If a named list
#'   of annotations is provided, uses the list names.
#'
#' @param sample_metadata Optional data frame with columns \code{marker},
#'   \code{sample_id}, and \code{replicate}. Metadata can either be derived from
#'   \code{bw_files} when filenames follow the expected naming convention, or supplied
#'   explicitly by the user for non-standard naming. One row per BigWig file is required;
#'   if a \code{bw_file} column is present it is used for filename-based matching,
#'   otherwise rows are matched by order. Default: \code{NULL} (auto-extract).
#'
#' @param bigwig_scale Character; which BigWig scaling to use. One of
#'   \code{"unscaled"}, \code{"scaled"}, or \code{"both"}. Default:
#'   \code{"unscaled"}.
#'
#' @param replicate_mode Character; which replicate type to include. One of
#'   \code{"all"}, \code{"pooled"}, or \code{"replicates"}. Default:
#'   \code{"all"}.
#'
#' @param scaling_info_file Optional path to a scaling info file (for example,
#'   \code{scalinginfo.txt}) containing at least \code{map_id} and \code{msr}
#'   columns. If provided, this file is used in preference to auto-discovery
#'   from \code{pipeline_output_path}. Default: \code{NULL}.
#'
#' @return An \code{EPK} object (S3 class) with slots:
#'   \itemize{
#'     \item \code{mse}: \code{MultiAssayExperiment} with one \code{SummarizedExperiment} per experiment
#'     \item \code{tables}: List containing \code{stats_summary}
#'     \item \code{enrichment_results}: List (empty by default; can be populated separately)
#'     \item \code{provenance}: Creation metadata (timestamp, session info)
#'   }
#'
#' @details
#' \strong{Annotation input flexibility:}
#' \itemize{
#'   \item \strong{GRanges}: Used directly. Must have metadata column \code{feature_id}
#'     or equivalent for row naming.
#'   \item \strong{BED path}: Read via \code{rtracklayer::import()}.
#'   \item \strong{Data frame}: Must have \code{chr}, \code{start}, \code{end} columns.
#'     \code{feature_id} column (or similar) is used for row names in the output.
#' }
#'
#' \strong{Marker extraction:}
#' Markers are detected from BigWig file names by looking for common patterns
#' (H3K4me3, H3K27ac, 5mC, etc.). Replicates are identified and handled. For
#' non-standard file names, provide \code{sample_metadata} with the expected columns
#' instead of relying on automatic parsing.
#'
#' \strong{Consistency checkpoint:}
#' If \code{stats_summary} is available, each BigWig file must follow
#' \code{(.*_(pooled|rep[0-9]+).genome.[unscaled|scaled].bw)} and its basename without
#' \code{.scaled/.unscaled.bw} must be present in \code{stats_summary$map_id}.
#' The function stops with an error if mismatches are found.
#'
#' \strong{BigWig filtering:}
#' Use \code{bigwig_scale} and \code{replicate_mode} to subset discovered/input
#' BigWig files before validation and signal extraction.
#'
#' \strong{Scaling info integration (path mode):}
#' When \code{pipeline_output_path} is provided, \code{create_epk()} attempts to
#' read \code{scalinginfo.txt} (or \code{scaling_info.txt}) from common
#' \code{reports/} locations and append \code{msr} to
#' \code{epk$tables$stats_summary} by \code{map_id} when available. If
#' \code{scaling_info_file} is provided, that file is used instead.
#'
#' @examples
#' \dontrun{
#' # Example 1: Explicit mode with GRanges annotation
#' data(toy_metadata, toy_genes)
#' toy_dir <- system.file("extdata", "toy_dataset", package = "EpigenicR")
#' bw_files <- list.files(
#'   file.path(toy_dir, "minute_output", "bigwig"),
#'   pattern = "\\.bw$",
#'   full.names = TRUE
#' )
#' stats_summary <- read.table(
#'   file.path(toy_dir, "minute_output", "reports", "stats_summary.txt"),
#'   header = TRUE,
#'   sep = "\t",
#'   stringsAsFactors = FALSE
#' )
#'
#' epk <- create_epk(
#'   bw_files = bw_files,
#'   annotations = toy_genes,
#'   stats_summary = stats_summary
#' )
#'
#' # Example 2: Path mode with BED annotation
#' epk <- create_epk(
#'   pipeline_output_path = "/path/to/project",
#'   annotations = "/path/to/annotations.bed"
#' )
#'
#' # Example 3: Multiple annotation sets
#' epk <- create_epk(
#'   bw_files = bw_files,
#'   annotations = list(
#'     genes = toy_genes,
#'     cpg_islands = cpg_islands_granges
#'   ),
#'   stats_summary = stats_summary
#' )
#'
#' # Example 4: Non-standard naming with explicit sample_metadata
#' metadata <- data.frame(
#'   marker = c("INPUT", "INPUT", "H3K4me3", "H3K4me3"),
#'   sample_id = c("Donor1_S1", "Donor1_S2", "Donor1_S1", "Donor1_S2"),
#'   replicate = c("pooled", "rep1", "pooled", "rep1")
#' )
#' epk <- create_epk(
#'   bw_files = bw_files,
#'   annotations = toy_genes,
#'   stats_summary = stats_summary,
#'   sample_metadata = metadata
#' )
#' }
#'
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @export
create_epk <- function(
  bw_files = NULL,
  annotations = NULL,
  stats_summary = NULL,
  pipeline_output_path = NULL,
  genome = "hg38",
  markers_to_exclude = c("INPUT"),
  experiment_names = NULL,
  sample_metadata = NULL,
  bigwig_scale = c("unscaled", "scaled", "both"),
  replicate_mode = c("all", "pooled", "replicates"),
  scaling_info_file = NULL,
  label_by = c("sample_id", "sample_id_batch")
) {

  # ---------- local helpers ----------
  .prepare_bw_metadata <- function(bw_files, sample_metadata = NULL) {
    bw_base <- basename(bw_files)

    if (is.null(sample_metadata)) {
      out <- .extract_marker_metadata(bw_files)
      out$bw_file <- bw_base
      out <- out[, c("bw_file", "marker", "sample_id", "replicate"), drop = FALSE]

      bad <- is.na(out$marker) | out$marker == "" |
        is.na(out$sample_id) | out$sample_id == "" |
        is.na(out$replicate) | out$replicate == ""

      if (any(bad)) {
        bad_files <- out$bw_file[bad]
        template <- data.frame(
          bw_file = head(bw_base, min(6, length(bw_base))),
          marker = "<marker>",
          sample_id = "<sample_id>",
          replicate = "<replicate: rep1/rep2/pooled>",
          stringsAsFactors = FALSE
        )

        message("Could not reliably parse marker/sample_id/replicate from some BigWig filenames.")
        message("Examples of problematic files:")
        message("  ", paste(head(bad_files, 10), collapse = "\n  "))
        message("Please provide 'sample_metadata' with columns: marker, sample_id, replicate (optional: bw_file).")
        message("Example template:")
        print(template)

        stop(
          "Automatic metadata extraction failed for one or more BigWig files. ",
          "Provide 'sample_metadata' explicitly."
        )
      }

      return(out)
    }

    required_cols <- c("marker", "sample_id", "replicate")
    if (!all(required_cols %in% names(sample_metadata))) {
      stop(
        "'sample_metadata' must have columns: ",
        paste(required_cols, collapse = ", ")
      )
    }

    sm <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)
    sm$marker <- trimws(as.character(sm$marker))
    sm$sample_id <- trimws(as.character(sm$sample_id))
    sm$replicate <- trimws(as.character(sm$replicate))

    if ("bw_file" %in% names(sm)) {
      sm$bw_file <- basename(trimws(as.character(sm$bw_file)))

      if (anyDuplicated(sm$bw_file)) {
        dup <- unique(sm$bw_file[duplicated(sm$bw_file)])
        stop(
          "'sample_metadata$bw_file' contains duplicates:\n  ",
          paste(dup, collapse = "\n  ")
        )
      }

      missing_meta <- setdiff(bw_base, sm$bw_file)
      extra_meta <- setdiff(sm$bw_file, bw_base)

      if (length(missing_meta) > 0) {
        stop(
          "'sample_metadata$bw_file' is missing entries for:\n  ",
          paste(missing_meta, collapse = "\n  ")
        )
      }
      if (length(extra_meta) > 0) {
        stop(
          "'sample_metadata$bw_file' has entries not present in 'bw_files':\n  ",
          paste(extra_meta, collapse = "\n  ")
        )
      }

      ord <- match(bw_base, sm$bw_file)
      out <- sm[ord, c("bw_file", "marker", "sample_id", "replicate"), drop = FALSE]
    } else {
      if (nrow(sm) != length(bw_files)) {
        stop(
          "'sample_metadata' must have one row per BigWig file when 'bw_file' is absent. ",
          "Expected ", length(bw_files), " rows, got ", nrow(sm), "."
        )
      }
      out <- data.frame(
        bw_file = bw_base,
        marker = sm$marker,
        sample_id = sm$sample_id,
        replicate = sm$replicate,
        stringsAsFactors = FALSE
      )
    }

    for (col in c("marker", "sample_id", "replicate")) {
      bad <- is.na(out[[col]]) | out[[col]] == ""
      if (any(bad)) {
        stop("Invalid 'sample_metadata': column '", col, "' contains NA/empty values.")
      }
    }

    out
  }

  .enrich_stats_summary_from_bw_metadata <- function(stats_summary, bw_metadata, replicate_mode = "all") {
    if (is.null(stats_summary) || !"map_id" %in% names(stats_summary)) {
      return(stats_summary)
    }

    bw_map_id <- sub("\\.(unscaled|scaled)\\.bw$", "", bw_metadata$bw_file)

    md <- data.frame(
      map_id = bw_map_id,
      marker = bw_metadata$marker,
      sample_id = bw_metadata$sample_id,
      replicate = bw_metadata$replicate,
      stringsAsFactors = FALSE
    )

    if (anyDuplicated(md$map_id)) {
      md <- md[!duplicated(md$map_id), , drop = FALSE]
      warning("Duplicate map_id values detected in bw metadata (e.g., scaled/unscaled pairs). Using first occurrence per map_id.")
    }

    # When replicate_mode is "pooled", bw files are named _pooled but
    # stats_summary tracks them as _rep1 (e.g. for cNUC product).
    # Normalise stats_summary map_ids to _pooled for matching.
    if (replicate_mode == "pooled") {
      stats_lookup_id <- sub("_rep[0-9]+", "_pooled", stats_summary$map_id)
    } else {
      stats_lookup_id <- stats_summary$map_id
    }

    ord <- match(stats_lookup_id, md$map_id)

    if (!"marker" %in% names(stats_summary)) {
      stats_summary$marker <- md$marker[ord]
    }
    if (!"sample_id" %in% names(stats_summary)) {
      stats_summary$sample_id <- md$sample_id[ord]
    }
    if (!"replicate" %in% names(stats_summary)) {
      stats_summary$replicate <- md$replicate[ord]
    }
    if (!"sample_id_rep" %in% names(stats_summary)) {
      stats_summary$sample_id_rep <- ifelse(
        is.na(stats_summary$sample_id) | is.na(stats_summary$replicate),
        NA_character_,
        paste(stats_summary$sample_id, stats_summary$replicate, sep = "_")
      )
    }

    stats_summary
  }

  .filter_bw_files <- function(bw_files, bigwig_scale, replicate_mode) {
    bw_base <- basename(bw_files)

    scale <- tolower(sub(".*\\.(scaled|unscaled)\\.bw$", "\\1", bw_base))
    scale[!grepl("\\.(scaled|unscaled)\\.bw$", bw_base, ignore.case = TRUE)] <- NA_character_

    rep_token <- tolower(sub(
      ".*_(pooled|rep[0-9]+)\\.[^.]+\\.(scaled|unscaled)\\.bw$",
      "\\1",
      bw_base,
      perl = TRUE
    ))
    rep_token[!grepl(
      "_(pooled|rep[0-9]+)\\.[^.]+\\.(scaled|unscaled)\\.bw$",
      bw_base,
      perl = TRUE,
      ignore.case = TRUE
    )] <- NA_character_

    keep <- rep(TRUE, length(bw_files))

    if (bigwig_scale != "both") {
      keep <- keep & !is.na(scale) & scale == bigwig_scale
    }

    if (replicate_mode == "pooled") {
      keep <- keep & !is.na(rep_token) & rep_token == "pooled"
    } else if (replicate_mode == "replicates") {
      keep <- keep & !is.na(rep_token) & grepl("^rep[0-9]+$", rep_token)
    }

    filtered <- bw_files[keep]
    message(
      "Selected ", length(filtered), " of ", length(bw_files),
      " BigWig files after filtering (bigwig_scale='", bigwig_scale,
      "', replicate_mode='", replicate_mode, "')."
    )

    if (length(filtered) == 0) {
      stop(
        "No BigWig files remain after filtering. Adjust 'bigwig_scale' and/or ",
        "'replicate_mode'."
      )
    }

    filtered
  }

  # ===== INPUT VALIDATION =====

  bigwig_scale <- match.arg(bigwig_scale)
  replicate_mode <- match.arg(replicate_mode)
  label_by <- match.arg(label_by)
  
  # Check that annotations is provided
  if (is.null(annotations)) {
    stop("'annotations' is required. Provide either a GRanges, BED file path, ",
         "data frame, or named list of these.")
  }

  # Determine input mode
  if (!is.null(pipeline_output_path) && !is.null(bw_files)) {
    warning(
      "'bw_files' is ignored because 'pipeline_output_path' is provided. ",
      "Using path-based mode."
    )
  }

  # ===== MODE 1: PIPELINE OUTPUT PATH =====
  if (!is.null(pipeline_output_path)) {
    message("Using path-based mode: discovering files from '", pipeline_output_path, "'")
    
    # Discover BigWig files
    bw_files <- .discover_bigwig_files(pipeline_output_path)
    if (length(bw_files) == 0) {
      stop("No BigWig files found in pipeline output path: ", pipeline_output_path)
    }
    message("Found ", length(bw_files), " BigWig files before filtering.")

    # Discover stats_summary
    stats_summary <- .discover_stats_summary(pipeline_output_path)
    if (is.null(stats_summary)) {
      warning(
        "No stats_summary found in pipeline output path. ",
        "Proceeding with minimal metadata; enrichment functions may be limited."
      )
    }
  } else {
    # ===== MODE 2: EXPLICIT INPUT =====
    # Validate explicit inputs
    if (is.null(bw_files)) {
      stop(
        "Either 'bw_files' (explicit mode) or 'pipeline_output_path' (path mode) ",
        "must be provided."
      )
    }

    if (!all(file.exists(bw_files))) {
      missing <- bw_files[!file.exists(bw_files)]
      stop("The following BigWig files do not exist:\n  ",
           paste(missing, collapse = "\n  "))
    }

    if (is.null(stats_summary)) {
      warning(
        "No 'stats_summary' provided. Proceeding with minimal metadata; ",
        "enrichment functions may be limited."
      )
    }
  }

  # Filter input files by scaling and replicate mode before validation.
  bw_files <- .filter_bw_files(
    bw_files = bw_files,
    bigwig_scale = bigwig_scale,
    replicate_mode = replicate_mode
  )

  if (!is.null(stats_summary)) {
    .validate_bw_files_in_stats_summary(
      bw_files = bw_files,
      stats_summary = stats_summary,
      replicate_mode = replicate_mode
    )
  }

  # ===== NORMALIZE ANNOTATIONS =====
  message("Normalizing annotations...")
  if (is.null(experiment_names)) {
    # Single annotation set
    if (!is.list(annotations) || inherits(annotations, c("GRanges", "data.frame"))) {
      experiment_names <- "primary_annotation"
      annotations <- list(primary_annotation = annotations)
    } else if (is.list(annotations) && !is.null(names(annotations))) {
      # Named list assume names are experiment names
      experiment_names <- names(annotations)
    } else {
      stop(
        "'annotations' must be a GRanges, BED path, data frame, or NAMED list. ",
        "Provide 'experiment_names' if you want to use an unnamed list."
      )
    }
  } else {
    # User provided experiment_names
    if (!is.list(annotations)) {
      annotations <- list(annotations)
      names(annotations) <- experiment_names
    } else if (!all(names(annotations) == experiment_names)) {
      names(annotations) <- experiment_names
    }
  }

  # Normalize each annotation to GRanges
  annotations <- lapply(
    annotations,
    .normalize_annotation_to_granges,
    genome = genome
  )

  # ===== EXTRACT MARKER NAMES FROM BIGWIG FILES =====
  message("Extracting marker information from BigWig filenames...")
  bw_metadata <- .prepare_bw_metadata(bw_files = bw_files, sample_metadata = sample_metadata)

  if (!is.null(sample_metadata)) {
    message("Using user-provided sample metadata.")
  }

  # Enrich stats_summary for downstream plotting helpers
  stats_summary <- .enrich_stats_summary_from_bw_metadata(
    stats_summary, bw_metadata, replicate_mode = replicate_mode
  )

  if (!is.null(stats_summary)) {
    stats_summary <- as.data.frame(stats_summary, stringsAsFactors = FALSE)

    if (all(c("raw_mapped", "raw_demultiplexed") %in% names(stats_summary))) {
      num <- suppressWarnings(as.numeric(stats_summary$raw_mapped))
      den <- suppressWarnings(as.numeric(stats_summary$raw_demultiplexed))
      stats_summary$frac_mapped <- ifelse(
        is.na(den) | den == 0,
        NA_real_,
        num / den
      )
    }

    if ("frac_mapq_filtered" %in% names(stats_summary)) {
      stats_summary$frac_mapq_filtered <- NULL
    }

    stats_summary <- .attach_msr_from_scaling_info(
      stats_summary = stats_summary,
      pipeline_output_path = pipeline_output_path,
      scaling_info_file = scaling_info_file
    )
  }

  # Determine unique markers to process
  unique_markers <- setdiff(unique(bw_metadata$marker), markers_to_exclude)
  unique_markers <- unique_markers[!is.na(unique_markers) & unique_markers != ""]
  if (length(unique_markers) == 0) {
    stop("No markers to process after exclusions. Check 'markers_to_exclude' and metadata.")
  }
  message("Processing markers: ", paste(unique_markers, collapse = ", "))

  # Enforce same sample set across markers for valid SummarizedExperiment assays
  marker_samples <- lapply(unique_markers, function(m) {
    idx <- !is.na(bw_metadata$marker) & bw_metadata$marker == m
    sort(unique(stats::na.omit(bw_metadata$sample_id[idx])))
  })
  names(marker_samples) <- unique_markers
  ref <- marker_samples[[1]]
  same_set <- vapply(marker_samples, function(x) identical(x, ref), logical(1))
  if (!all(same_set)) {
    details <- vapply(
      names(marker_samples),
      function(m) paste0(m, " [n=", length(marker_samples[[m]]), "]: ", paste(marker_samples[[m]], collapse = ", ")),
      character(1)
    )
    stop(
      "Markers do not share the same sample set, so assays would have different ncol.\n",
      "Provide balanced files/metadata so each marker has identical sample_id values.\n",
      paste(details, collapse = "\n")
    )
  }

  # ===== BUILD EXPERIMENTS =====
  message("Building SummarizedExperiment objects...")
  experiments <- list()

  for (exp_name in experiment_names) {
    message("  Processing experiment: ", exp_name)
    annotation_gr <- annotations[[exp_name]]

    # Build assay list (one per marker)
    assay_list <- list()

    for (marker in unique_markers) {
      message("    Marker: ", marker)
      marker_idx <- !is.na(bw_metadata$marker) & bw_metadata$marker == marker
      marker_bw_files <- bw_files[marker_idx]
      marker_labels <- if (label_by == "sample_id_batch") {
        paste(bw_metadata$sample_id[marker_idx],
              bw_metadata$batch[marker_idx], sep = "_")
      } else {
        bw_metadata$sample_id[marker_idx]
      }

      if (length(marker_bw_files) == 0) {
        stop("No BigWig files available for marker '", marker, "' after filtering.")
      }

      # Use wigglescout to extract signal at regions
      bw_gr <- .extract_bw_signal(
        bwfiles = marker_bw_files,
        loci = annotation_gr,
        labels = marker_labels
      )

      # Convert to matrix
      m <- as.matrix(S4Vectors::mcols(bw_gr)[, , drop = FALSE])
      storage.mode(m) <- "numeric"

      assay_list[[marker]] <- m
    }

    # Create SummarizedExperiment
    coldata <- S4Vectors::DataFrame(
      sample_id = colnames(assay_list[[1]]),
      row.names = colnames(assay_list[[1]])
    )

    se <- SummarizedExperiment::SummarizedExperiment(
      assays = assay_list,
      rowRanges = annotation_gr,
      colData = coldata
    )

    # Set row names based on annotation feature IDs
    rownames(se) <- .get_feature_ids(annotation_gr)

    experiments[[exp_name]] <- se
  }

  # ===== BUILD MULTIASSAYEXPERIMENT =====
  message("Building MultiAssayExperiment...")
  coldata_mse <- S4Vectors::DataFrame(
    sample_id = colnames(experiments[[1]]),
    row.names = colnames(experiments[[1]])
  )

  mse <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experiments,
    colData = coldata_mse
  )

  # ===== BUILD EPK OBJECT =====
  message("Creating EPK object...")

  tables_list <- list(stats_summary = stats_summary)

  epk <- structure(
    list(
      mse = mse,
      tables = tables_list,
      enrichment_results = list(
        chromatin_states = NULL,
        enrichment_profile = NULL
      ),
      provenance = list(
        created = Sys.time(),
        session = sessionInfo()
      )
    ),
    class = "EPK"
  )

  message("EPK object successfully created!")
  return(epk)
}


# ===== HELPER FUNCTIONS =====

#' Discover BigWig files in pipeline output directory
#' @keywords internal
.discover_bigwig_files <- function(path) {
  # Search for .bw files in minute-chip style directories first
  candidates <- c(
    file.path(path, "minute_output", "bigwig"),
    file.path(path, "bigwig"),
    file.path(path, "bw"),
    file.path(path, "BigWig"),
    path  # or in root
  )

  for (dir in candidates) {
    if (dir.exists(dir)) {
      bw_files <- list.files(dir, pattern = "\\.bw$", full.names = TRUE)
      if (length(bw_files) > 0) {
        return(bw_files)
      }
    }
  }

  return(character(0))
}

#' Discover stats_summary file in pipeline output directory
#' @keywords internal
.discover_stats_summary <- function(path) {
  candidates <- c(
    file.path(path, "minute_output", "reports", "stats_summary.txt"),
    file.path(path, "minute_output", "reports", "stats_summary.csv"),
    file.path(path, "minute_output", "reports", "stats_summary", "stats_summary.txt"),
    file.path(path, "minute_output", "reports", "stats_summary", "stats_summary.csv"),
    file.path(path, "reports", "stats_summary.txt"),
    file.path(path, "reports", "stats_summary.csv"),
    file.path(path, "reports", "stats_summary", "stats_summary.txt"),
    file.path(path, "reports", "stats_summary", "stats_summary.csv"),
    file.path(path, "stats_summary.txt"),
    file.path(path, "stats_summary.csv"),
    file.path(path, "qc", "stats_summary.txt"),
    file.path(path, "qc", "stats_summary.csv")
  )

  for (file in candidates) {
    if (file.exists(file)) {
      sep <- if (grepl("\\.csv$", file)) "," else "\t"
      return(read.csv(file, sep = sep, stringsAsFactors = FALSE))
    }
  }

  return(NULL)
}

#' Discover and attach msr from scalinginfo file
#' @keywords internal
.attach_msr_from_scaling_info <- function(stats_summary, pipeline_output_path, scaling_info_file = NULL) {
  if (is.null(stats_summary)) {
    return(stats_summary)
  }

  if (!is.null(scaling_info_file)) {
    if (!is.character(scaling_info_file) || length(scaling_info_file) != 1) {
      stop("'scaling_info_file' must be a single file path string.")
    }
    if (!file.exists(scaling_info_file)) {
      message("No scaling info file found at '", scaling_info_file, "'. Continuing without msr.")
      return(stats_summary)
    }
    scaling_file <- scaling_info_file
  } else {
    if (is.null(pipeline_output_path)) {
      return(stats_summary)
    }

    candidates <- c(
      file.path(pipeline_output_path, "minute_output", "reports", "scalinginfo.txt"),
      file.path(pipeline_output_path, "reports", "scalinginfo.txt"),
      file.path(pipeline_output_path, "minute_output", "reports", "scaling_info.txt"),
      file.path(pipeline_output_path, "reports", "scaling_info.txt")
    )

    scaling_file <- candidates[file.exists(candidates)][1]
    if (is.na(scaling_file) || length(scaling_file) == 0) {
      message("No scaling info file found in pipeline reports directories. Continuing without msr.")
      return(stats_summary)
    }
  }

  scaling_info <- tryCatch(
    utils::read.csv(
      scaling_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      warning("Could not read scaling info file: ", scaling_file, " (", conditionMessage(e), ").")
      return(NULL)
    }
  )

  if (is.null(scaling_info)) {
    return(stats_summary)
  }

  if (!"map_id" %in% names(scaling_info) && ncol(scaling_info) >= 1) {
    names(scaling_info)[1] <- "map_id"
  }

  if (!all(c("map_id", "msr") %in% names(scaling_info))) {
    warning("Scaling info found but missing required columns 'map_id' and/or 'msr': ", scaling_file)
    return(stats_summary)
  }

  ord <- match(stats_summary$map_id, scaling_info$map_id)
  stats_summary$msr <- scaling_info$msr[ord]

  stats_summary
}

#' Validate that BigWig files match stats_summary$map_id
#' @keywords internal
.validate_bw_files_in_stats_summary <- function(bw_files, stats_summary,
                                                replicate_mode = "all") {
  if (!"map_id" %in% names(stats_summary)) {
    stop("'stats_summary' must contain a 'map_id' column.")
  }

  bw_base <- basename(bw_files)
  bw_pattern <- "_(pooled|rep[0-9]+)\\.[^.]+\\.(unscaled|scaled)\\.bw$"
  invalid_bw <- bw_base[!grepl(bw_pattern, bw_base)]
  if (length(invalid_bw) > 0) {
    stop(
      "Invalid BigWig filename format. Expected pattern: ",
      "(.*_(pooled|rep[0-9]+).genome.(unscaled|scaled).bw)\n",
      "Invalid file(s):\n  ",
      paste(invalid_bw, collapse = "\n  ")
    )
  }

  map_id_pattern <- "_(pooled|rep[0-9]+)\\.[^.]+$"
  invalid_map_id <- stats_summary$map_id[!grepl(map_id_pattern, stats_summary$map_id)]
  if (length(invalid_map_id) > 0) {
    stop(
      "Invalid map_id format in stats_summary. Expected pattern: ",
      "(.*_(pooled|rep[0-9]+).genome)\n",
      "Invalid map_id value(s):\n  ",
      paste(unique(invalid_map_id), collapse = "\n  ")
    )
  }

  bw_map_id <- sub("\\.(unscaled|scaled)\\.bw$", "", bw_base)

  # When replicate_mode is "pooled", bw files are named _pooled but
  # stats_summary tracks them as _rep1 (e.g. for cNUC product).
  # Normalise stats_summary map_ids to _pooled for matching.
  if (replicate_mode == "pooled") {
    stats_summary_map_id <- sub("_rep[0-9]+", "_pooled", stats_summary$map_id)
  } else {
    stats_summary_map_id <- stats_summary$map_id
  }

  missing_in_stats <- setdiff(bw_map_id, stats_summary_map_id)

  if (length(missing_in_stats) > 0) {
    stop(
      "The following bw_files are not listed in stats_summary$map_id:\n  ",
      paste(missing_in_stats, collapse = "\n  ")
    )
  }

  invisible(TRUE)
}

#' Normalize annotation input to GRanges
#' @keywords internal
.normalize_annotation_to_granges <- function(annotation, genome) {
  # If already GRanges, return as-is
  if (inherits(annotation, "GRanges")) {
    return(annotation)
  }

  # If file path, read it
  if (is.character(annotation) && length(annotation) == 1 && file.exists(annotation)) {
    message("    Reading annotation from file: ", annotation)
    return(rtracklayer::import(annotation, format = "BED"))
  }

  # If data frame, convert to GRanges
  if (is.data.frame(annotation)) {
    message("    Converting data frame to GRanges")
    required_cols <- c("chr", "start", "end")
    if (!all(required_cols %in% names(annotation))) {
      stop(
        "Data frame annotation must have columns: ",
        paste(required_cols, collapse = ", ")
      )
    }

    gr <- GenomicRanges::GRanges(
      seqnames = annotation$chr,
      ranges = IRanges::IRanges(annotation$start, annotation$end),
      genome = genome
    )

    # Add metadata columns if present
    mcols(gr) <- annotation[, setdiff(names(annotation), required_cols)]

    return(gr)
  }

  stop(
    "Annotation must be a GRanges object, BED file path, or data frame. ",
    "Got class: ", class(annotation)
  )
}

#' Extract BigWig signal at loci
#' @keywords internal
.extract_bw_signal <- function(bwfiles, loci, labels) {
  # Wrapper around wigglescout::bw_loci
  # For now, assume wigglescout is available

  if (!requireNamespace("wigglescout", quietly = TRUE)) {
    stop(
      "Package 'wigglescout' is required but not installed. ",
      "Install it from Bioconductor: BiocManager::install('wigglescout')"
    )
  }

  bw_gr <- wigglescout::bw_loci(
    bwfiles = bwfiles,
    loci = loci,
    labels = labels
  )

  return(bw_gr)
}

#' Extract marker names and metadata from BigWig filenames
#' @keywords internal
.extract_marker_metadata <- function(bw_files) {
  # Pattern: Project_Batch_Marker_Rerun_Sample_Rep.Genome.Scaling.bw
  # e.g.: Proj1_A1_H3K4me3_1_SAMPLE-0008_pooled.hg38.unscaled.bw

  metadata <- data.frame(
    bw_file   = basename(bw_files),
    marker    = NA_character_,
    batch     = NA_character_,
    sample_id = NA_character_,
    replicate = NA_character_,
    stringsAsFactors = FALSE
  )

  # Define known markers (can be extended)
  known_markers <- c(
    "H3K4me3", "H3K9me3", "H3K27me3", "H3K27ac",
    "5mC", "CXXC", "INPUT"
  )

  for (i in seq_along(bw_files)) {
    filename <- basename(bw_files[i])

    # Primary parser for minute-style naming:
    # <project>_<batch>_<marker>_<rerun>_<sample_id>_<replicate>.<genome>.<scaled|unscaled>.bw
    m <- regexec(
      "^([^_]+)_([^_]+)_([^_]+)_([^_]+)_(.+)_(pooled|rep[0-9]+)\\.[^.]+\\.(scaled|unscaled)\\.bw$",
      filename,
      ignore.case = TRUE,
      perl = TRUE
    )
    g <- regmatches(filename, m)[[1]]

    if (length(g) > 0) {
      metadata$marker[i]    <- g[4]
      metadata$batch[i]     <- g[3]
      metadata$sample_id[i] <- g[6]
      metadata$replicate[i] <- tolower(g[7])
      next
    }

    # Extract marker
    for (marker in known_markers) {
      if (grepl(marker, filename, fixed = TRUE)) {
        metadata$marker[i] <- marker
        break
      }
    }

    # Fallback: infer marker as the third underscore-delimited token.
    if (is.na(metadata$marker[i]) && grepl("^[^_]+_[^_]+_[^_]+_", filename)) {
      metadata$marker[i] <- sub("^[^_]+_[^_]+_([^_]+)_.*$", "\\1", filename)
    }

    # Extract sample_id (simple heuristic: look for SAMPLE-XXXX pattern)
    sample_match <- regmatches(
      filename,
      gregexpr("SAMPLE-[0-9]+", filename)
    )
    if (length(sample_match[[1]]) > 0) {
      metadata$sample_id[i] <- sample_match[[1]][1]
    } else {
      # Fallback: use first part of filename
      metadata$sample_id[i] <- sub("_.*", "", filename)
    }

    # Extract replicate (rep1, rep2, pooled, etc.)
    rep_match <- regmatches(
      filename,
      gregexpr("(rep[0-9]+|pooled)", filename, ignore.case = TRUE)
    )
    if (length(rep_match[[1]]) > 0) {
      metadata$replicate[i] <- tolower(rep_match[[1]][1])
    } else {
      metadata$replicate[i] <- "unknown"
    }
  }

  if (any(is.na(metadata$marker))) {
    warning(
      "Could not extract marker from some BigWig filenames. ",
      "Check file naming convention."
    )
  }

  return(metadata)
}

#' Get feature IDs from GRanges metadata
#' @keywords internal
.get_feature_ids <- function(gr) {
  # Try common metadata column names
  candidates <- c(
    "feature_id", "cpg_id", "peak_id", "region_id",
    "gene_name", "gene_id", "id", "name"
  )

  for (col in candidates) {
    if (col %in% names(mcols(gr))) {
      return(mcols(gr)[[col]])
    }
  }

  # Fallback: use default names
  return(paste0("region_", seq_along(gr)))
}
