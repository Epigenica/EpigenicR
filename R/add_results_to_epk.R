#' Add enrichment and ChromHMM results to an EPK object
#'
#' Reads pre-computed enrichment profile and chromatin state distribution CSV
#' files from a results directory and stores them in the corresponding slots of
#' an existing \code{EPK} object.
#'
#' @details
#' Two directory layouts are supported and can be mixed within the same
#' \code{results_path}:
#'
#' \describe{
#'   \item{Nested layout}{Marker subdirectories sit inside an annotation
#'     directory, e.g. \code{protein_coding/H3K27ac/H3K27ac_profile_start_data.csv}.
#'     The subdirectory name is used as the marker key.}
#'   \item{Flat layout}{CSV files sit directly inside an annotation directory,
#'     e.g. \code{CpG_islands/methylation_profile_start_data.csv}.
#'     The marker key is inferred from the filename prefix up to the first
#'     \code{_profile} or \code{_chromatin_state_dist} token.}
#' }
#'
#' Files are routed to slots based on their filename:
#' \itemize{
#'   \item \code{*_chromatin_state_dist.csv} →
#'     \code{epk$enrichment_results$chromatin_states[[annotation]][[marker]]}
#'   \item \code{*_profile_*_data.csv} →
#'     \code{epk$enrichment_results$enrichment_profile[[annotation]][[marker]]}
#' }
#'
#' @param epk An \code{EPK} object.
#' @param results_path Character; path to the top-level results directory whose
#'   immediate subdirectories represent annotation sets (e.g.
#'   \code{"protein_coding"}, \code{"CpG_islands"}).
#'
#' @return The updated \code{EPK} object with
#'   \code{enrichment_results$chromatin_states} and
#'   \code{enrichment_results$enrichment_profile} populated.
#'
#' @examples
#' \dontrun{
#' epk <- add_results_to_epk(epk, results_path = "path/to/results")
#' str(epk$enrichment_results, max.level = 2)
#' }
#'
#' @export
add_results_to_epk <- function(epk, results_path) {
  stopifnot(inherits(epk, "EPK"))
  if (!dir.exists(results_path)) {
    stop("'results_path' does not exist: ", results_path)
  }

  .read_csv_safe <- function(path) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # Infer marker from filename prefix when there is no marker subdirectory,
  # e.g. "methylation_profile_start_data.csv" -> "methylation"
  .marker_from_filename <- function(fname) {
    sub(
      "_(profile|chromatin_state_dist).*",
      "",
      tools::file_path_sans_ext(fname)
    )
  }

  .slot_file <- function(epk, csv_path, annotation, marker) {
    fname <- basename(csv_path)
    data  <- .read_csv_safe(csv_path)

    if (grepl("_chromatin_state_dist\\.csv$", fname)) {
      epk$enrichment_results$chromatin_states[[annotation]][[marker]] <- data
    } else if (grepl("_profile_.*_data\\.csv$", fname)) {
      epk$enrichment_results$enrichment_profile[[annotation]][[marker]] <- data
    } else {
      message("Skipping unrecognised file: ", fname)
    }
    epk
  }

  annotation_dirs <- list.dirs(results_path, full.names = TRUE, recursive = FALSE)

  if (length(annotation_dirs) == 0) {
    warning("No annotation subdirectories found in: ", results_path)
    return(epk)
  }

  for (ann_dir in annotation_dirs) {
    annotation <- basename(ann_dir)

    # Nested layout: immediate subdirs are marker directories
    marker_dirs <- list.dirs(ann_dir, full.names = TRUE, recursive = FALSE)
    for (mk_dir in marker_dirs) {
      marker    <- basename(mk_dir)
      csv_files <- list.files(mk_dir, pattern = "\\.csv$", full.names = TRUE)
      for (f in csv_files) {
        epk <- .slot_file(epk, f, annotation, marker)
      }
    }

    # Flat layout: CSV files directly in annotation dir (no marker subdir)
    csv_files <- list.files(
      ann_dir,
      pattern    = "\\.csv$",
      full.names = TRUE,
      recursive  = FALSE
    )
    for (f in csv_files) {
      marker <- .marker_from_filename(basename(f))
      epk    <- .slot_file(epk, f, annotation, marker)
    }
  }

  epk
}
