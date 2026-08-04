#' Add enrichment and ChromHMM results to an EPK object
#'
#' Reads pre-computed enrichment profile and chromatin state distribution CSV
#' files from a results directory and stores them in the corresponding slots of
#' an existing \code{EPK} object.
#'
#' @details
#' Profile results and chromatin state results live in separate root directories
#' and are loaded with two separate calls:
#'
#' \describe{
#'   \item{Profiles}{\code{output/profile/<loci_name>/<marker>/} — keyed as
#'     \code{enrichment_profile[[loci_name]][[marker]]}.}
#'   \item{Chromatin states}{\code{output/chromhmm/<chromhmm_ref>/<marker>/} — keyed as
#'     \code{chromatin_states[[chromhmm_ref]][[marker]]}, where \code{chromhmm_ref}
#'     is the ChromHMM reference BED filename without the \code{.bed} extension
#'     (e.g. \code{"E107_15_coreMarks_hg38lift_mnemonics"}).}
#' }
#'
#' Within each root, two sub-layouts are supported:
#' \describe{
#'   \item{Nested}{\code{<key>/<marker>/file.csv} — marker inferred from subdirectory name.}
#'   \item{Flat}{\code{<key>/file.csv} — marker inferred from filename prefix.}
#' }
#'
#' Files are routed to slots based on their filename:
#' \itemize{
#'   \item \code{*_chromatin_state_dist.csv} →
#'     \code{epk$enrichment_results$chromatin_states[[chromhmm_ref]][[marker]]}
#'   \item \code{*_profile_*_data.csv} →
#'     \code{epk$enrichment_results$enrichment_profile[[loci_name]][[marker]]}
#' }
#'
#' @param epk An \code{EPK} object.
#' @param results_path Character; path to a top-level results directory whose
#'   immediate subdirectories are used as the slot key (loci name for profiles,
#'   ChromHMM reference name for chromatin states).
#'
#' @return The updated \code{EPK} object with the relevant
#'   \code{enrichment_results} slots populated.
#'
#' @examples
#' \dontrun{
#' # Load profiles (keyed by loci name)
#' epk <- add_results_to_epk(epk, results_path = "output/profile")
#'
#' # Load chromatin states (keyed by ChromHMM reference name)
#' epk <- add_results_to_epk(epk, results_path = "output/chromhmm")
#'
#' # Results are now accessible:
#' # epk$enrichment_results$enrichment_profile$protein_coding$H3K4me3
#' # epk$enrichment_results$chromatin_states$E107_15_coreMarks_hg38lift_mnemonics$H3K4me3
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
