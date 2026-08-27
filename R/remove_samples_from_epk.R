#' Remove samples from every slot of an EPK object
#'
#' Drops the given samples from \code{epk$mse} (every experiment, \code{colData},
#' and \code{sampleMap}), from any data frame with a \code{sample} column or
#' square named matrix nested inside \code{epk$enrichment_results} and
#' \code{epk$derived}, and from \code{epk$tables$stats_summary} /
#' \code{epk$tables$sample_metadata}.
#'
#' @param epk An \code{EPK} object.
#' @param samples Character vector of sample identifiers to remove, given
#'   exactly as they appear in \code{colnames(epk$mse)} (e.g.
#'   \code{"A.c_rep1"}) — not a bare \code{sample_id}. This is deliberate:
#'   a biological sample can span multiple replicate columns, and the exact
#'   column-naming convention (replicate suffix, batch suffix, R's
#'   \code{make.names()} mangling of digit-leading IDs) varies enough across
#'   projects that guessing which columns belong to a sample is unsafe. Use
#'   \code{grep("^<sample_id>", colnames(epk$mse), value = TRUE)} to find the
#'   exact strings first.
#'
#' @return The updated \code{EPK} object with matching samples removed.
#'
#' @details
#' \strong{What gets checked where:}
#' \itemize{
#'   \item \code{epk$mse}: subset via \code{MultiAssayExperiment}'s own
#'     primary-based \code{[} method, so every experiment, \code{colData}, and
#'     \code{sampleMap} stay consistent. Samples not found in
#'     \code{colData(epk$mse)} are skipped with a warning.
#'   \item \code{epk$enrichment_results} / \code{epk$derived}: walked
#'     recursively. Any data frame with a column literally named
#'     \code{sample} has matching rows dropped; any matrix with row/column
#'     names has the matching row \emph{and} column dropped (for sample x
#'     sample correlation matrices). Anything else (including \code{NULL}
#'     slots) is left untouched.
#'   \item \code{epk$tables$stats_summary}: rows filtered by the
#'     \code{sample_id_rep} column only — \emph{not} \code{sample_id}. If a
#'     row's \code{sample_id_rep} uses yet another naming variant than what
#'     you pass to \code{samples}, it will not be removed; filter
#'     \code{epk$tables$stats_summary} by \code{sample_id} directly for that
#'     case.
#'   \item \code{epk$tables$sample_metadata}: keyed by base \code{sample_id}
#'     (no replicate suffix). A row is only dropped once none of its
#'     replicate columns remain in \code{epk$mse} after removal.
#' }
#'
#' \strong{Known gotcha — digit-leading sample IDs:} purely numeric
#' \code{sample_id} values (e.g. \code{"7156"}) are not syntactically valid R
#' names, so \code{make.names()} prepends an \code{"X"} wherever it gets
#' applied (e.g. some \code{epk$mse} column-naming paths) — but not
#' everywhere (e.g. \code{enrichment_profile}'s \code{sample} column, written
#' directly from CSV output, keeps the original unprefixed value). The same
#' sample can therefore appear as both \code{"X7156_rep1"} and
#' \code{"7156_rep1"} in different slots of the same object. This function
#' does not try to reconcile that automatically — check each slot's actual
#' values first (e.g. \code{grep("7156", ..., value = TRUE)}) and call this
#' function once per naming variant if needed.
#'
#' @examples
#' \dontrun{
#' epk <- remove_samples_from_epk(epk, samples = c("A.c_rep1", "A.c_rep2"))
#' }
#'
#' @importFrom SummarizedExperiment colData
#' @export
remove_samples_from_epk <- function(epk, samples) {
  stopifnot(inherits(epk, "EPK"))
  samples <- as.character(samples)

  # ---- epk$mse: MultiAssayExperiment's own primary-based subsetting keeps
  # every experiment, colData, and sampleMap consistent ----
  primaries <- rownames(SummarizedExperiment::colData(epk$mse))
  missing <- setdiff(samples, primaries)
  if (length(missing) > 0) {
    warning(
      "Not found in epk$mse colData, skipping there: ",
      paste(missing, collapse = ", ")
    )
  }
  keep <- setdiff(primaries, samples)
  epk$mse <- epk$mse[, keep, ]

  # ---- recursively drop matching rows from any data.frame with a $sample
  # column, and matching row+col from any square named matrix ----
  .scrub <- function(x) {
    if (is.data.frame(x)) {
      if ("sample" %in% names(x)) {
        return(x[!(x$sample %in% samples), , drop = FALSE])
      }
      return(x)
    }
    if (is.matrix(x)) {
      rn <- rownames(x)
      cn <- colnames(x)
      if (!is.null(rn) && !is.null(cn) && any(rn %in% samples)) {
        return(x[setdiff(rn, samples), setdiff(cn, samples), drop = FALSE])
      }
      return(x)
    }
    if (is.list(x)) {
      return(lapply(x, .scrub))
    }
    x
  }

  epk$enrichment_results <- .scrub(epk$enrichment_results)
  epk$derived <- .scrub(epk$derived)

  # ---- epk$tables$stats_summary: per-replicate rows ----
  if (!is.null(epk$tables$stats_summary) &&
      "sample_id_rep" %in% names(epk$tables$stats_summary)) {
    epk$tables$stats_summary <- epk$tables$stats_summary[
      !(epk$tables$stats_summary$sample_id_rep %in% samples), ,
      drop = FALSE
    ]
  }

  # ---- epk$tables$sample_metadata: keyed by base sample_id; only drop once
  # no remaining epk$mse columns still reference it ----
  if (!is.null(epk$tables$sample_metadata) &&
      "sample_id" %in% names(epk$tables$sample_metadata)) {
    remaining_base_ids <- sub("_(pooled|rep[0-9]+)$", "", keep)
    epk$tables$sample_metadata <- epk$tables$sample_metadata[
      epk$tables$sample_metadata$sample_id %in% remaining_base_ids, ,
      drop = FALSE
    ]
  }

  epk$provenance$updated <- Sys.time()
  epk
}
