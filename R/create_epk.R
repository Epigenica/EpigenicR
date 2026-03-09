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
#'   \code{(.*_rep[0-9].genome)} so BigWig files can be checked against
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
#' (H3K4me3, H3K27ac, 5mC, etc.). Replicates are identified and handled.
#'
#' \strong{Consistency checkpoint:}
#' If \code{stats_summary} is available, each BigWig file must follow
#' \code{(.*_rep[0-9].genome.[unscaled|scaled].bw)} and its basename without
#' \code{.scaled/.unscaled.bw} must be present in \code{stats_summary$map_id}.
#' The function stops with an error if mismatches are found.
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
  experiment_names = NULL
) {

  # ===== INPUT VALIDATION =====
  
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
    message("Found ", length(bw_files), " BigWig files.")

    # Discover stats_summary
    stats_summary <- .discover_stats_summary(pipeline_output_path)
    if (is.null(stats_summary)) {
      warning(
        "No stats_summary found in pipeline output path. ",
        "Proceeding with minimal metadata; enrichment functions may be limited."
      )
    } else {
      .validate_bw_files_in_stats_summary(bw_files = bw_files, stats_summary = stats_summary)
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
    } else {
      .validate_bw_files_in_stats_summary(bw_files = bw_files, stats_summary = stats_summary)
    }
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
  bw_metadata <- .extract_marker_metadata(bw_files)

  # Determine unique markers to process
  unique_markers <- setdiff(unique(bw_metadata$marker), markers_to_exclude)
  if (length(unique_markers) == 0) {
    stop("No markers to process after exclusions. Check 'markers_to_exclude'.")
  }
  message("Processing markers: ", paste(unique_markers, collapse = ", "))

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
      marker_bw_files <- bw_files[bw_metadata$marker == marker]

      # Use wigglescout to extract signal at regions
      bw_gr <- .extract_bw_signal(
        bwfiles = marker_bw_files,
        loci = annotation_gr,
        labels = bw_metadata$sample_id[bw_metadata$marker == marker]
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

#' Validate that BigWig files match stats_summary$map_id
#' @keywords internal
.validate_bw_files_in_stats_summary <- function(bw_files, stats_summary) {
  if (!"map_id" %in% names(stats_summary)) {
    stop("'stats_summary' must contain a 'map_id' column.")
  }

  bw_base <- basename(bw_files)
  bw_pattern <- "_rep[0-9]+\\.[^.]+\\.(unscaled|scaled)\\.bw$"
  invalid_bw <- bw_base[!grepl(bw_pattern, bw_base)]
  if (length(invalid_bw) > 0) {
    stop(
      "Invalid BigWig filename format. Expected pattern: ",
      "(.*_rep[0-9].genome.[unscaled|scaled].bw)\n",
      "Invalid file(s):\n  ",
      paste(invalid_bw, collapse = "\n  ")
    )
  }

  map_id_pattern <- "_rep[0-9]+\\.[^.]+$"
  invalid_map_id <- stats_summary$map_id[!grepl(map_id_pattern, stats_summary$map_id)]
  if (length(invalid_map_id) > 0) {
    stop(
      "Invalid map_id format in stats_summary. Expected pattern: ",
      "(.*_rep[0-9].genome)\n",
      "Invalid map_id value(s):\n  ",
      paste(unique(invalid_map_id), collapse = "\n  ")
    )
  }

  bw_map_id <- sub("\\.(unscaled|scaled)\\.bw$", "", bw_base)
  missing_in_stats <- setdiff(bw_map_id, stats_summary$map_id)

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
    bw_file = basename(bw_files),
    marker = NA_character_,
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

    # Extract marker
    for (marker in known_markers) {
      if (grepl(marker, filename, fixed = TRUE)) {
        metadata$marker[i] <- marker
        break
      }
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
  candidates <- c("feature_id", "gene_name", "gene_id", "id", "name")

  for (col in candidates) {
    if (col %in% names(mcols(gr))) {
      return(mcols(gr)[[col]])
    }
  }

  # Fallback: use default names
  return(paste0("region_", seq_along(gr)))
}
