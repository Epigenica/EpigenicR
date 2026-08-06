#' Plot BigWig enrichment profile at genomic loci
#'
#' Runs \code{wigglescout::plot_bw_profile} for a set of pre-resolved BigWig
#' files, saves the resulting plot and underlying data to \code{output_dir},
#' and returns \code{invisible(NULL)}.
#'
#' This function is intentionally agnostic to marker type — the caller chooses
#' \code{mode} and \code{loci_label} appropriate for their context (e.g.
#' \code{mode = "start"} for histone TSS profiles, \code{mode = "center"} for
#' methylation at CpG islands).
#'
#' @param allfiles Character vector of BigWig file paths.
#' @param allfiles_name Character vector of sample labels corresponding to
#'   \code{allfiles}.
#' @param loci A \code{GRanges} object specifying the genomic features.
#' @param mk Character; marker name used in output filenames and the x-axis
#'   label (e.g. \code{"H3K4me3"}, \code{"5mC"}).
#' @param output_dir Character; directory where PNG and CSV files are written.
#'   Created recursively if absent.
#' @param mode Character; anchoring mode passed to
#'   \code{wigglescout::plot_bw_profile}. One of \code{"start"} (default),
#'   \code{"stretch"}, \code{"end"}, or \code{"center"}.
#' @param loci_label Character; short description of the loci type shown in
#'   the x-axis label, e.g. \code{"Protein coding"} or \code{"CpG islands"}.
#'   Default: \code{"genomic features"}.
#'
#' @return Invisibly returns \code{NULL}. Writes to \code{output_dir}:
#'   \itemize{
#'     \item PNG: \code{<mk>_profile_<mode>.png}
#'     \item CSV: \code{<mk>_profile_<mode>_data.csv}
#'   }
#'
#' @seealso \code{\link{run_chromhmm_histone}}, \code{\link{run_chromhmm_methylation}},
#'   \code{\link{dispatch_jobs}}
#'
#' @examples
#' \dontrun{
#' # Histone — TSS anchored at start
#' run_bw_profile(
#'   allfiles      = c("/bw/H3K4me3_S1.bw", "/bw/INPUT_S1.bw"),
#'   allfiles_name = c("H3K4me3_S1_rep1",   "INPUT_S1_rep1"),
#'   loci          = tss_granges,
#'   mk            = "H3K4me3",
#'   output_dir    = "results/protein_coding/H3K4me3",
#'   mode          = "start",
#'   loci_label    = "Protein coding"
#' )
#'
#' # Methylation — centered on CpG islands
#' run_bw_profile(
#'   allfiles      = c("/bw/5mC_S1.bw", "/bw/CXXC_S1.bw"),
#'   allfiles_name = c("5mC_S1_rep1",   "CXXC_S1_rep1"),
#'   loci          = cpg_granges,
#'   mk            = "5mC",
#'   output_dir    = "results/CpG_islands/5mC",
#'   mode          = "center",
#'   loci_label    = "CpG islands"
#' )
#' }
#'
#' @export
run_bw_profile <- function(allfiles, allfiles_name, loci, mk, output_dir,
                           mode = c("start", "stretch", "end", "center"),
                           loci_label = "genomic features") {
  mode <- match.arg(mode)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  p_enrich <- wigglescout::plot_bw_profile(
    allfiles,
    loci   = loci,
    mode   = mode,
    labels = allfiles_name
  )

  p_enrich <- p_enrich +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::labs(
      x = paste0(mk, " ", loci_label, " (", length(loci), ") loci")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 1, byrow = TRUE)) +
    ggplot2::theme(
      legend.position = "right",
      legend.title    = ggplot2::element_blank()
    )

  profile_png <- file.path(output_dir, paste0(mk, "_profile_", mode, ".png"))
  ggplot2::ggsave(profile_png, plot = p_enrich, width = 12, height = 8)

  profile_csv <- file.path(
    output_dir,
    paste0(mk, "_profile_", mode, "_data.csv")
  )
  write.table(
    p_enrich$data,
    profile_csv,
    row.names = FALSE,
    quote     = FALSE,
    col.names = TRUE,
    sep       = ","
  )

  invisible(NULL)
}


#' Compute ChromHMM state distribution for a histone modification marker
#'
#' Filters BigWig metadata for a histone marker (plus INPUT controls), computes
#' per-chromatin-state mean RPGC via
#' \code{wigglescout::plot_bw_loci_summary_heatmap}, and saves the summary
#' table as a CSV. Profile generation is intentionally excluded — call
#' \code{\link{run_bw_profile}} separately if a profile is also needed.
#'
#' @param bw_df Data frame of BigWig metadata with columns \code{marker},
#'   \code{sample_id}, \code{replicate}, \code{batch}, and \code{bw_file}.
#' @param bigwig_dir Character; directory containing the BigWig files listed
#'   in \code{bw_df$bw_file}.
#' @param mk Character; marker name to process (e.g. \code{"H3K4me3"},
#'   \code{"H3K27ac"}).
#' @param output_dir Character; directory for output files (created if absent).
#' @param chromHmm_path Character; path to the ChromHMM annotation directory.
#' @param chromHMM_annotation Character; ChromHMM annotation filename within
#'   \code{chromHmm_path} (e.g.
#'   \code{"E107_15_coreMarks_hg38lift_mnemonics.bed"}).
#' @param product Character; \code{"cNUC"} selects unscaled BigWigs; all other
#'   values select scaled BigWigs. INPUT controls are always included.
#' @param replicate_type Character; \code{"replicate"} (default, non-pooled)
#'   or \code{"pooled"}.
#'
#' @return Invisibly returns \code{NULL}. Writes to \code{output_dir}:
#'   \itemize{
#'     \item CSV: \code{<mk>_chromatin_state_dist.csv} with columns
#'       \code{Chromatin_State}, \code{sample_id_rep}, \code{mean_rpgc_val},
#'       \code{mean_rpgc_text}
#'     \item Sentinel: \code{.done} on successful completion
#'   }
#'
#' @details
#' Chromatin state enrichment is genome-wide against the ChromHMM reference
#' BED — no user-supplied loci are involved.
#'
#' \strong{Output directory convention:} use
#' \code{output/chromhmm/<chromhmm_ref>/<marker>} where \code{chromhmm_ref}
#' is the ChromHMM annotation filename without \code{.bed}
#' (e.g. \code{gsub(".bed$", "", chromHMM_annotation)}). This ensures
#' \code{add_results_to_epk()} keys results by ChromHMM reference rather than
#' by loci name, allowing multiple references to coexist in one EPK.
#'
#' Designed for parallel execution via \code{\link{dispatch_jobs}}.
#' The \code{.done} sentinel allows downstream cache-checking to skip completed
#' markers on re-runs.
#'
#' For methylation markers (5mC / CXXC), use
#' \code{\link{run_chromhmm_methylation}} instead, which includes CXXC in the
#' file selection.
#'
#' @seealso \code{\link{run_bw_profile}}, \code{\link{run_chromhmm_methylation}},
#'   \code{\link{dispatch_jobs}}
#'
#' @examples
#' \dontrun{
#' chromhmm_ref <- gsub(".bed$", "", "E107_15_coreMarks_hg38lift_mnemonics.bed")
#' run_chromhmm_histone(
#'   bw_df               = bw_df,
#'   bigwig_dir          = "/path/to/bigwigs",
#'   mk                  = "H3K4me3",
#'   output_dir          = file.path("output/chromhmm", chromhmm_ref, "H3K4me3"),
#'   chromHmm_path       = "data/chromHmm_annotations",
#'   chromHMM_annotation = "E107_15_coreMarks_hg38lift_mnemonics.bed",
#'   product             = "chromatin"
#' )
#' }
#'
#' @export
run_chromhmm_histone <- function(bw_df, bigwig_dir, mk, output_dir,
                                 chromHmm_path, chromHMM_annotation,
                                 product,
                                 replicate_type = c("replicate", "pooled")) {
  replicate_type <- match.arg(replicate_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (replicate_type == "pooled") {
    bw_df_subset <- dplyr::filter(
      bw_df,
      replicate == "pooled",
      marker == mk | tolower(marker) == "input"
    )
  } else {
    bw_df_subset <- dplyr::filter(
      bw_df,
      replicate != "pooled",
      marker == mk | tolower(marker) == "input"
    )
  }

  if (product != "cNUC") {
    bw_df_subset <- dplyr::filter(
      bw_df_subset,
      grepl("\\.scaled\\.", bw_file, perl = TRUE) | tolower(marker) == "input"
    )
  } else {
    bw_df_subset <- dplyr::filter(
      bw_df_subset,
      grepl("\\.unscaled", bw_file) | tolower(marker) == "input"
    )
  }

  allfiles <- file.path(bigwig_dir, bw_df_subset$bw_file)

  if (length(allfiles) == 0) {
    message(sprintf("[chromHMM:%s] no BigWig files found in bw_df - skipping.", mk))
    file.create(file.path(output_dir, ".done"))
    return(invisible(NULL))
  }

  if (length(unique(bw_df_subset$batch)) > 1) {
    allfiles_name <- paste0(
      bw_df_subset$marker, "_", bw_df_subset$sample_id,
      "_", bw_df_subset$replicate, "_", bw_df_subset$batch
    )
  } else {
    allfiles_name <- paste0(
      bw_df_subset$marker, "_", bw_df_subset$sample_id,
      "_", bw_df_subset$replicate
    )
  }

  p_heatmap <- wigglescout::plot_bw_loci_summary_heatmap(
    allfiles,
    file.path(chromHmm_path, chromHMM_annotation),
    labels     = allfiles_name,
    remove_top = 0.01
  )

  tmp_df <- p_heatmap@data
  colnames(tmp_df) <- c("Chromatin_State", "sample_id_rep", "mean_rpgc_val", "mean_rpgc_text")

  write.table(
    tmp_df,
    file.path(output_dir, paste0(mk, "_chromatin_state_dist.csv")),
    sep       = ",",
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  file.create(file.path(output_dir, ".done"))
  invisible(NULL)
}


#' Compute ChromHMM state distribution for methylation markers
#'
#' Filters BigWig metadata for methylation markers — including the paired CXXC
#' control — computes per-chromatin-state mean RPGC via
#' \code{wigglescout::plot_bw_loci_summary_heatmap}, and saves the summary
#' table as a CSV. Profile generation is intentionally excluded — call
#' \code{\link{run_bw_profile}} separately if a profile is also needed.
#'
#' Only call this function when the project contains methylation markers (5mC).
#' For projects with histone marks only, use \code{\link{run_chromhmm_histone}}.
#'
#' @param bw_df Data frame of BigWig metadata with columns \code{marker},
#'   \code{sample_id}, \code{replicate}, \code{batch}, and \code{bw_file}.
#' @param bigwig_dir Character; directory containing the BigWig files listed
#'   in \code{bw_df$bw_file}.
#' @param mk Character; primary methylation marker (typically \code{"5mC"}).
#'   CXXC files are included automatically.
#' @param output_dir Character; directory for output files (created if absent).
#' @param chromHmm_path Character; path to the ChromHMM annotation directory.
#' @param chromHMM_annotation Character; ChromHMM annotation filename within
#'   \code{chromHmm_path}.
#' @param product Character; \code{"cNUC"} selects unscaled BigWigs; all other
#'   values select scaled BigWigs. INPUT controls are always included.
#' @param replicate_type Character; \code{"replicate"} (default, non-pooled)
#'   or \code{"pooled"}.
#'
#' @return Invisibly returns \code{NULL}. Writes to \code{output_dir}:
#'   \itemize{
#'     \item CSV: \code{<mk>_chromatin_state_dist.csv} with columns
#'       \code{Chromatin_State}, \code{sample_id_rep}, \code{mean_rpgc_val},
#'       \code{mean_rpgc_text}
#'     \item Sentinel: \code{.done} on successful completion
#'   }
#'
#' @details
#' CXXC files are pulled from \code{bw_df} alongside \code{mk} so that the
#' ChromHMM enrichment reflects both 5mC and its unmodified counterpart in the
#' same run. Chromatin state enrichment is genome-wide against the ChromHMM
#' reference BED — no user-supplied loci are involved.
#'
#' \strong{Output directory convention:} use
#' \code{output/chromhmm/<chromhmm_ref>/<marker>} — see
#' \code{\link{run_chromhmm_histone}} for details.
#'
#' Designed for parallel execution via \code{\link{dispatch_jobs}}.
#'
#' @seealso \code{\link{run_bw_profile}}, \code{\link{run_chromhmm_histone}},
#'   \code{\link{dispatch_jobs}}
#'
#' @examples
#' \dontrun{
#' # Only call when 5mC is present in the project
#' chromhmm_ref <- gsub(".bed$", "", "E107_15_coreMarks_hg38lift_mnemonics.bed")
#' run_chromhmm_methylation(
#'   bw_df               = bw_df,
#'   bigwig_dir          = "/path/to/bigwigs",
#'   mk                  = "5mC",
#'   output_dir          = file.path("output/chromhmm", chromhmm_ref, "5mC"),
#'   chromHmm_path       = "data/chromHmm_annotations",
#'   chromHMM_annotation = "E107_15_coreMarks_hg38lift_mnemonics.bed",
#'   product             = "chromatin"
#' )
#' }
#'
#' @export
run_chromhmm_methylation <- function(bw_df, bigwig_dir, mk, output_dir,
                                     chromHmm_path, chromHMM_annotation,
                                     product,
                                     replicate_type = c("replicate", "pooled")) {
  replicate_type <- match.arg(replicate_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (replicate_type == "pooled") {
    bw_df_subset <- dplyr::filter(
      bw_df,
      replicate == "pooled",
      marker == mk | marker == "CXXC" | tolower(marker) == "input"
    )
  } else {
    bw_df_subset <- dplyr::filter(
      bw_df,
      replicate != "pooled",
      marker == mk | marker == "CXXC" | tolower(marker) == "input"
    )
  }

  if (product != "cNUC") {
    bw_df_subset <- dplyr::filter(
      bw_df_subset,
      grepl("\\.scaled\\.", bw_file, perl = TRUE) | tolower(marker) == "input"
    )
  } else {
    bw_df_subset <- dplyr::filter(
      bw_df_subset,
      grepl("\\.unscaled", bw_file) | tolower(marker) == "input"
    )
  }

  allfiles <- file.path(bigwig_dir, bw_df_subset$bw_file)

  if (length(allfiles) == 0) {
    message(sprintf("[chromHMM:%s] no BigWig files found in bw_df - skipping.", mk))
    file.create(file.path(output_dir, ".done"))
    return(invisible(NULL))
  }

  if (length(unique(bw_df_subset$batch)) > 1) {
    allfiles_name <- paste0(
      bw_df_subset$marker, "_", bw_df_subset$sample_id,
      "_", bw_df_subset$replicate, "_", bw_df_subset$batch
    )
  } else {
    allfiles_name <- paste0(
      bw_df_subset$marker, "_", bw_df_subset$sample_id,
      "_", bw_df_subset$replicate
    )
  }

  p_heatmap <- wigglescout::plot_bw_loci_summary_heatmap(
    allfiles,
    file.path(chromHmm_path, chromHMM_annotation),
    labels     = allfiles_name,
    remove_top = 0.01
  )

  tmp_df <- p_heatmap@data
  colnames(tmp_df) <- c("Chromatin_State", "sample_id_rep", "mean_rpgc_val", "mean_rpgc_text")

  write.table(
    tmp_df,
    file.path(output_dir, paste0(mk, "_chromatin_state_dist.csv")),
    sep       = ",",
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  file.create(file.path(output_dir, ".done"))
  invisible(NULL)
}


#' Dispatch and orchestrate parallel ChromHMM jobs
#'
#' Manages parallel execution of ChromHMM state distribution jobs using
#' background R processes via \code{callr}. Maintains a worker pool, dispatches
#' pending jobs, and monitors completion.
#'
#' @param jobs List of job specifications, each a named list with:
#'   \itemize{
#'     \item \code{fn}: Function to execute
#'       (\code{\link{run_chromhmm_histone}} or
#'       \code{\link{run_chromhmm_methylation}}).
#'     \item \code{args}: Named list of arguments to pass to \code{fn}.
#'     \item \code{mk}: Character; marker name used for logging.
#'   }
#' @param n_workers Integer; maximum number of parallel workers.
#'
#' @return Invisibly returns \code{NULL} after all jobs complete.
#'
#' @details
#' \enumerate{
#'   \item Maintains a pool of up to \code{n_workers} background processes.
#'   \item Polls running processes; upon completion checks exit status and logs.
#'   \item Dispatches pending jobs to idle workers.
#'   \item Warns on non-zero exit; messages on success.
#'   \item Loops until all jobs are complete.
#' }
#'
#' @seealso \code{\link{run_bw_profile}}, \code{\link{run_chromhmm_histone}}, \code{\link{run_chromhmm_methylation}}
#'
#' @examples
#' \dontrun{
#' histone_markers <- c("H3K4me3", "H3K27ac", "H3K9me3")
#' chromhmm_annotation <- "E107_15_coreMarks_hg38lift_mnemonics.bed"
#' chromhmm_ref        <- gsub(".bed$", "", chromhmm_annotation)
#'
#' jobs <- lapply(histone_markers, function(mk) {
#'   list(
#'     fn   = run_chromhmm_histone,
#'     mk   = mk,
#'     args = list(
#'       bw_df               = bw_df,
#'       bigwig_dir          = "/path/to/bw",
#'       mk                  = mk,
#'       output_dir          = file.path("output/chromhmm", chromhmm_ref, mk),
#'       chromHmm_path       = "data/chromHmm_annotations",
#'       chromHMM_annotation = chromhmm_annotation,
#'       product             = "chromatin"
#'     )
#'   )
#' })
#'
#' dispatch_jobs(jobs, n_workers = 4)
#'
#' # Methylation separately, only if 5mC is in the project
#' if ("5mC" %in% unique(bw_df$marker)) {
#'   run_chromhmm_methylation(bw_df = bw_df, mk = "5mC", ...)
#' }
#' }
#'
#' @export
dispatch_jobs <- function(jobs, n_workers) {
  running <- list()
  pending <- jobs

  repeat {
    done <- vapply(running, function(j) !j$proc$is_alive(), logical(1))

    for (i in which(done)) {
      j   <- running[[i]]
      err <- j$proc$read_all_error()

      if (j$proc$get_exit_status() != 0) {
        warning(sprintf(
          "[chromHMM] worker '%s' FAILED (exit %d):\n%s",
          j$mk, j$proc$get_exit_status(), err
        ))
      } else {
        if (nzchar(trimws(err))) {
          message(sprintf("[chromHMM] worker '%s' stderr:\n%s", j$mk, err))
        }
        message(sprintf("[chromHMM] completed: %s", j$mk))
      }
    }

    running <- running[!done]

    while (length(running) < n_workers && length(pending) > 0) {
      job          <- pending[[1]]
      pending[[1]] <- NULL

      proc <- callr::r_bg(
        func = function(fn, fn_args) {
          library(GenomeInfoDb,  quietly = TRUE)  # nolint: undesirable_function_linter
          library(GenomicRanges, quietly = TRUE)  # nolint: undesirable_function_linter
          library(wigglescout,   quietly = TRUE)  # nolint: undesirable_function_linter
          library(ggplot2,       quietly = TRUE)  # nolint: undesirable_function_linter
          library(dplyr,         quietly = TRUE)  # nolint: undesirable_function_linter
          library(stringr,       quietly = TRUE)  # nolint: undesirable_function_linter
          do.call(fn, fn_args)
        },
        args = list(fn = job$fn, fn_args = job$args)
      )
      running <- c(running, list(list(proc = proc, mk = job$mk)))

      message(sprintf(
        "[chromHMM] dispatched: %s (%d/%d running)",
        job$mk, length(running), n_workers
      ))
    }

    if (length(running) == 0 && length(pending) == 0) break

    Sys.sleep(1)
  }

  invisible(NULL)
}
