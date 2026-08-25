#' Add new marker assays to an existing EPK object
#'
#' Extends an existing \code{EPK} object by adding one or more new marker assays
#' to every \code{SummarizedExperiment} already present in \code{epk$mse}.
#' Genomic regions are taken directly from the EPK's stored \code{rowRanges},
#' so annotations do not need to be re-supplied.
#'
#' @param epk Existing \code{EPK} object created by \code{\link{create_epk}}.
#'
#' @param bw_files Character vector of BigWig file paths (explicit mode).
#'   Ignored if \code{pipeline_output_path} is provided.
#'
#' @param pipeline_output_path Character path to pipeline output directory
#'   (path mode). If provided, \code{bw_files} is ignored.
#'
#' @param stats_summary Optional data frame with QC statistics. Falls back to
#'   \code{epk$tables$stats_summary} when \code{NULL} and no
#'   \code{pipeline_output_path} is given.
#'
#' @param sample_metadata Optional data frame with columns \code{marker},
#'   \code{sample_id}, and \code{replicate}. Use for non-standard filenames.
#'   See \code{\link{create_epk}} for details.
#'
#' @param genome Character; genome version (default \code{"hg38"}).
#'
#' @param markers_to_exclude Character vector of marker names to skip
#'   (default \code{c("INPUT")}). Matching is exact and case-sensitive. To
#'   include \code{INPUT} as its own assay, set
#'   \code{markers_to_exclude = character(0)}.
#'
#' @param markers_to_include Character vector of marker names to process.
#'   When \code{NULL} (default) all detected markers minus exclusions are used.
#'
#' @param bigwig_scale Character; scaling tier to load. One of
#'   \code{"unscaled"}, \code{"scaled"}, or \code{"both"}
#'   (default \code{"unscaled"}).
#'
#' @param replicate_mode Character; replicates to include. One of
#'   \code{"all"} (default), \code{"pooled"}, or \code{"replicates"}.
#'
#' @param scaling_info_file Optional path to a scaling info file containing
#'   \code{map_id} and \code{msr} columns.
#'
#' @param label_by Character; sample labelling. One of \code{"sample_id"}
#'   (default) or \code{"sample_id_batch"}.
#'
#' @param overwrite Logical; if \code{FALSE} (default) the function stops when
#'   a new marker assay name already exists in any experiment. Set \code{TRUE}
#'   to replace existing assays.
#'
#' @return Updated \code{EPK} object with new assays added to every experiment
#'   in \code{epk$mse}.
#'
#' @details
#' Unlike \code{\link{add_features_to_epk}}, which adds new \emph{experiments}
#' (annotation sets / row sets) to the \code{MultiAssayExperiment}, this
#' function adds new \emph{assays} (markers) to every existing
#' \code{SummarizedExperiment} while keeping the row and column structure
#' intact.
#'
#' Samples in the new BigWig files must match those already in the EPK.
#' An error is raised for any mismatch; use \code{sample_metadata} to provide
#' explicit sample IDs when auto-detection is unreliable.
#'
#' @examples
#' \dontrun{
#' # EPK already contains 5mC, H3K27me3, H3K4me3, H3K9me3 assays.
#' # Add H3K27ac from a separate set of BigWig files:
#' epk <- add_marker_to_epk(
#'   epk              = epk,
#'   bw_files         = bw_files_h3k27ac,
#'   stats_summary    = stats_summary,
#'   markers_to_exclude = c("Input"),
#'   bigwig_scale     = "unscaled",
#'   replicate_mode   = "replicates"
#' )
#'
#' # Verify the new assay is present in every experiment
#' SummarizedExperiment::assayNames(epk$mse[["genes"]])
#' }
#'
#' @importFrom SummarizedExperiment rowRanges assay assayNames assay<-
#' @importFrom MultiAssayExperiment experiments MultiAssayExperiment
#' @importFrom SummarizedExperiment colData
#' @export
add_marker_to_epk <- function(
  epk,
  bw_files             = NULL,
  pipeline_output_path = NULL,
  stats_summary        = NULL,
  sample_metadata      = NULL,
  genome               = "hg38",
  markers_to_exclude   = c("INPUT"),
  markers_to_include   = NULL,
  bigwig_scale         = c("unscaled", "scaled", "both"),
  replicate_mode       = c("all", "pooled", "replicates"),
  scaling_info_file    = NULL,
  label_by             = c("sample_id", "sample_id_batch"),
  overwrite            = FALSE
) {
  bigwig_scale   <- match.arg(bigwig_scale)
  replicate_mode <- match.arg(replicate_mode)
  label_by       <- match.arg(label_by)

  if (is.null(epk) || !is.list(epk) || is.null(epk$mse)) {
    stop("'epk' must be a valid EPK object containing an 'mse' slot.")
  }

  if (is.null(pipeline_output_path) && is.null(bw_files)) {
    stop("Provide either 'pipeline_output_path' (path mode) or 'bw_files' (explicit mode).")
  }

  existing_exps <- MultiAssayExperiment::experiments(epk$mse)
  if (length(existing_exps) == 0) {
    stop("EPK contains no experiments. Use create_epk() to build from scratch.")
  }

  # Inherit stats_summary from EPK when not explicitly provided
  if (is.null(pipeline_output_path) && is.null(stats_summary) &&
      !is.null(epk$tables) && !is.null(epk$tables$stats_summary)) {
    stats_summary <- epk$tables$stats_summary
  }

  # Re-use the rowRanges already stored in the EPK as annotations
  annotations <- lapply(existing_exps, SummarizedExperiment::rowRanges)

  message("Adding new marker(s) to ", length(existing_exps), " experiment(s): ",
          paste(names(existing_exps), collapse = ", "))

  # Build a temporary EPK over the same regions with the new bw files
  epk_new <- create_epk(
    bw_files             = bw_files,
    annotations          = annotations,
    stats_summary        = stats_summary,
    sample_metadata      = sample_metadata,
    pipeline_output_path = pipeline_output_path,
    genome               = genome,
    markers_to_exclude   = markers_to_exclude,
    markers_to_include   = markers_to_include,
    bigwig_scale         = bigwig_scale,
    replicate_mode       = replicate_mode,
    scaling_info_file    = scaling_info_file,
    label_by             = label_by
  )

  new_exps         <- MultiAssayExperiment::experiments(epk_new$mse)
  existing_samples <- colnames(existing_exps[[1]])

  for (nm in names(existing_exps)) {
    se_existing <- existing_exps[[nm]]
    se_new      <- new_exps[[nm]]
    new_samples <- colnames(se_new)

    # Validate sample set
    if (!setequal(existing_samples, new_samples)) {
      missing_in_new <- setdiff(existing_samples, new_samples)
      extra_in_new   <- setdiff(new_samples, existing_samples)
      msg <- paste0("Sample mismatch in experiment '", nm, "'.")
      if (length(missing_in_new) > 0)
        msg <- paste0(msg, " Missing in new data: ", paste(missing_in_new, collapse = ", "), ".")
      if (length(extra_in_new) > 0)
        msg <- paste0(msg, " Extra in new data: ", paste(extra_in_new, collapse = ", "), ".")
      stop(msg)
    }

    # Align column order to the existing SE
    if (!identical(existing_samples, new_samples)) {
      se_new <- se_new[, existing_samples]
    }

    # Check for assay name conflicts
    new_assay_nms      <- SummarizedExperiment::assayNames(se_new)
    existing_assay_nms <- SummarizedExperiment::assayNames(se_existing)
    conflicts          <- intersect(new_assay_nms, existing_assay_nms)

    if (length(conflicts) > 0 && !isTRUE(overwrite)) {
      stop(
        "Assay(s) already exist in experiment '", nm, "': ",
        paste(conflicts, collapse = ", "),
        ". Set overwrite = TRUE to replace them."
      )
    }

    for (assay_nm in new_assay_nms) {
      SummarizedExperiment::assay(se_existing, assay_nm) <-
        SummarizedExperiment::assay(se_new, assay_nm)
    }

    existing_exps[[nm]] <- se_existing
  }

  epk$mse <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = existing_exps,
    colData     = SummarizedExperiment::colData(epk$mse)
  )
  epk$provenance$updated <- Sys.time()

  message("Done. New assay(s): ", paste(
    SummarizedExperiment::assayNames(epk$mse[[1]]), collapse = ", "
  ))

  return(epk)
}
