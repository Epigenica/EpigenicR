#' Add new feature experiments to an existing EPK object
#'
#' Convenience wrapper to extend an existing \\code{EPK} object with one or more
#' new feature experiments in \\code{mse}. Internally, this function builds the new
#' experiments using \\code{create_epk()} and merges them into the existing object.
#'
#' @param epk Existing \\code{EPK} object.
#' @param annotations Either a \\code{GRanges} object, path to a BED file, data frame,
#'   or a named list of such objects/paths for multiple feature sets.
#' @param pipeline_output_path Character path to pipeline output directory (path mode).
#' @param bw_files Character vector of BigWig file paths (explicit mode).
#' @param stats_summary Optional data frame with QC statistics. If not provided in
#'   explicit mode and available in \\code{epk$tables$stats_summary}, that table is used.
#' @param genome Character; genome version passed to \\code{create_epk()}.
#' @param markers_to_exclude Character vector of marker names to skip.
#' @param experiment_names Optional character vector naming added experiments.
#' @param overwrite Logical; if \\code{FALSE} (default), stop when experiment names
#'   already exist in \\code{epk$mse}. If \\code{TRUE}, replace existing experiments.
#'
#' @return Updated \\code{EPK} object with added/replaced experiments in \\code{mse}.
#'
#' @examples
#' \\dontrun{
#' toy_dir <- system.file("extdata", "toy_dataset", package = "EpigenicR")
#' data(toy_genes)
#' epk <- create_epk(pipeline_output_path = toy_dir, annotations = toy_genes)
#'
#' epk <- add_features_to_epk(
#'   epk = epk,
#'   pipeline_output_path = toy_dir,
#'   annotations = list(enhancer = "data/enhancer.bed")
#' )
#' }
#' @export
add_features_to_epk <- function(
  epk,
  annotations,
  pipeline_output_path = NULL,
  bw_files = NULL,
  stats_summary = NULL,
  genome = "hg38",
  markers_to_exclude = c("INPUT"),
  experiment_names = NULL,
  overwrite = FALSE
) {
  if (is.null(epk) || !is.list(epk) || is.null(epk$mse)) {
    stop("'epk' must be a valid EPK object containing an 'mse' slot.")
  }

  if (is.null(annotations)) {
    stop("'annotations' is required. Provide one or more feature annotations to add.")
  }

  if (is.null(pipeline_output_path) && is.null(bw_files)) {
    stop("Provide either 'pipeline_output_path' (path mode) or 'bw_files' (explicit mode).")
  }

  if (is.null(pipeline_output_path) && is.null(stats_summary) &&
      !is.null(epk$tables) && !is.null(epk$tables$stats_summary)) {
    stats_summary <- epk$tables$stats_summary
  }

  epk_new <- create_epk(
    bw_files = bw_files,
    annotations = annotations,
    stats_summary = stats_summary,
    pipeline_output_path = pipeline_output_path,
    genome = genome,
    markers_to_exclude = markers_to_exclude,
    experiment_names = experiment_names
  )

  existing_exps <- MultiAssayExperiment::experiments(epk$mse)
  new_exps <- MultiAssayExperiment::experiments(epk_new$mse)

  duplicate_names <- intersect(names(existing_exps), names(new_exps))
  if (length(duplicate_names) > 0 && !isTRUE(overwrite)) {
    stop(
      "Experiment name(s) already exist in epk$mse: ",
      paste(duplicate_names, collapse = ", "),
      ". Set overwrite = TRUE to replace them."
    )
  }

  if (length(existing_exps) > 0) {
    existing_samples <- colnames(existing_exps[[1]])

    for (nm in names(new_exps)) {
      new_samples <- colnames(new_exps[[nm]])

      if (!setequal(existing_samples, new_samples)) {
        missing_in_new <- setdiff(existing_samples, new_samples)
        extra_in_new <- setdiff(new_samples, existing_samples)

        msg <- paste0(
          "Sample mismatch while adding experiment '", nm, "'."
        )
        if (length(missing_in_new) > 0) {
          msg <- paste0(
            msg,
            " Missing in new data: ", paste(missing_in_new, collapse = ", "), "."
          )
        }
        if (length(extra_in_new) > 0) {
          msg <- paste0(
            msg,
            " Extra in new data: ", paste(extra_in_new, collapse = ", "), "."
          )
        }
        stop(msg)
      }

      new_exps[[nm]] <- new_exps[[nm]][, existing_samples]
    }
  }

  for (nm in names(new_exps)) {
    existing_exps[[nm]] <- new_exps[[nm]]
  }

  MultiAssayExperiment::experiments(epk$mse) <- existing_exps
  epk$provenance$updated <- Sys.time()

  return(epk)
}
