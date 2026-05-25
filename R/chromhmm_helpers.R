#' Plot BigWig enrichment profile at genomic loci
#'
#' Runs \code{wigglescout::plot_bw_profile} for a set of pre-resolved BigWig
#' files, saves the resulting plot and underlying data to \code{output_dir},
#' and returns \code{invisible(NULL)}. This is the shared profile engine called
#' by \code{\link{run_chromhmm_histone}} and \code{\link{run_chromhmm_methylation}},
#' but can also be called directly when file lists and labels are already built.
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
#' @seealso \code{\link{run_chromhmm_histone}}, \code{\link{run_chromhmm_methylation}}
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


#' Run BigWig enrichment profile and ChromHMM state distribution for histone markers
#'
#' Filters BigWig metadata for a histone marker (plus INPUT controls), resolves
#' file paths and sample labels, then:
#' \enumerate{
#'   \item Delegates profile plotting to \code{\link{run_bw_profile}}.
#'   \item Computes per-chromatin-state mean RPGC via
#'     \code{wigglescout::plot_bw_loci_summary_heatmap} and saves the summary
#'     table as a CSV.
#' }
#'
#' @param bw_df Data frame of BigWig metadata with columns \code{marker},
#'   \code{sample_id}, \code{replicate}, \code{batch}, and \code{bw_file}.
#' @param bigwig_dir Character; directory containing the BigWig files listed
#'   in \code{bw_df$bw_file}.
#' @param mk Character; marker name to process (e.g. \code{"H3K4me3"},
#'   \code{"H3K27ac"}).
#' @param loci A \code{GRanges} object of genomic features (e.g. protein-coding
#'   TSS regions).
#' @param output_dir Character; directory for output files (created if absent).
#' @param chromHmm_path Character; path to the ChromHMM annotation directory.
#' @param chromHMM_annotation Character; ChromHMM annotation filename within
#'   \code{chromHmm_path} (e.g.
#'   \code{"E107_15_coreMarks_hg38lift_mnemonics.bed"}).
#' @param product Character; product type. \code{"cNUC"} selects unscaled
#'   BigWigs; all other values select scaled BigWigs. INPUT controls are always
#'   included regardless of scaling.
#' @param replicate_type Character; which replicates to include.
#'   \code{"replicate"} (default) uses non-pooled replicates;
#'   \code{"pooled"} uses pooled files only.
#' @param plot_bw_profile_mode Character; anchoring mode passed to
#'   \code{wigglescout::plot_bw_profile}. One of \code{"start"} (default),
#'   \code{"stretch"}, \code{"end"}, or \code{"center"}.
#'
#' @return Invisibly returns \code{NULL}. Writes to \code{output_dir}:
#'   \itemize{
#'     \item PNG + CSV: enrichment profile via \code{\link{run_bw_profile}}
#'     \item CSV: chromatin state distribution
#'       (\code{<mk>_chromatin_state_dist.csv})
#'     \item Sentinel: \code{.done} on successful completion
#'   }
#'
#' @details
#' Designed for parallel execution via \code{\link{dispatch_chromhmm_jobs}}.
#' The \code{.done} sentinel allows \code{\link{run_chromhmm_enrichment}} to
#' skip completed markers on re-runs.
#'
#' @seealso \code{\link{run_bw_profile}}, \code{\link{run_chromhmm_methylation}},
#'   \code{\link{run_chromhmm_enrichment}}
#'
#' @examples
#' \dontrun{
#' bw_df <- data.frame(
#'   marker    = c("H3K4me3", "H3K4me3", "INPUT"),
#'   sample_id = c("Sample_A", "Sample_B", "Sample_A"),
#'   replicate = c("rep1", "rep1", "rep1"),
#'   batch     = c("B1", "B1", "B1"),
#'   bw_file   = c("file1.bw", "file2.bw", "file3.bw")
#' )
#' loci <- GenomicRanges::GRanges("chr1:1-1000:+")
#'
#' run_chromhmm_histone(
#'   bw_df               = bw_df,
#'   bigwig_dir          = "/path/to/bigwigs",
#'   mk                  = "H3K4me3",
#'   loci                = loci,
#'   output_dir          = "/path/to/output",
#'   chromHmm_path       = "/path/to/chromhmm",
#'   chromHMM_annotation = "E107_15_coreMarks_hg38lift_mnemonics.bed",
#'   product             = "chromatin"
#' )
#' }
#'
#' @export
run_chromhmm_histone <- function(bw_df, bigwig_dir, mk, loci, output_dir,
                                 chromHmm_path, chromHMM_annotation,
                                 product,
                                 replicate_type       = c("replicate", "pooled"),
                                 plot_bw_profile_mode = c("start", "stretch", "end", "center")) {
  replicate_type       <- match.arg(replicate_type)
  plot_bw_profile_mode <- match.arg(plot_bw_profile_mode)
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

  # --- enrichment profile ---
  run_bw_profile(
    allfiles      = allfiles,
    allfiles_name = allfiles_name,
    loci          = loci,
    mk            = mk,
    output_dir    = output_dir,
    mode          = plot_bw_profile_mode,
    loci_label    = "Protein coding"
  )

  # --- ChromHMM state distribution (data only, no boxplot) ---
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


#' Run BigWig enrichment profile and ChromHMM state distribution for methylation markers
#'
#' Filters BigWig metadata for methylation markers — including the paired CXXC
#' control — resolves file paths and sample labels, then:
#' \enumerate{
#'   \item Delegates profile plotting to \code{\link{run_bw_profile}} with
#'     \code{mode = "center"}.
#'   \item Computes per-chromatin-state mean RPGC via
#'     \code{wigglescout::plot_bw_loci_summary_heatmap} and saves the summary
#'     table as a CSV.
#' }
#'
#' @param bw_df Data frame of BigWig metadata with columns \code{marker},
#'   \code{sample_id}, \code{replicate}, \code{batch}, and \code{bw_file}.
#' @param bigwig_dir Character; directory containing the BigWig files listed
#'   in \code{bw_df$bw_file}.
#' @param mk Character; primary methylation marker (typically \code{"5mC"}).
#'   CXXC files are included automatically.
#' @param loci A \code{GRanges} object of genomic features (e.g. CpG islands).
#' @param output_dir Character; directory for output files (created if absent).
#' @param chromHmm_path Character; path to the ChromHMM annotation directory.
#' @param chromHMM_annotation Character; ChromHMM annotation filename within
#'   \code{chromHmm_path}.
#' @param product Character; product type. \code{"cNUC"} selects unscaled
#'   BigWigs; all other values select scaled BigWigs. INPUT controls are always
#'   included.
#' @param replicate_type Character; \code{"replicate"} (default, non-pooled) or
#'   \code{"pooled"}.
#'
#' @return Invisibly returns \code{NULL}. Writes to \code{output_dir}:
#'   \itemize{
#'     \item PNG + CSV: enrichment profile via \code{\link{run_bw_profile}}
#'       (always \code{mode = "center"})
#'     \item CSV: chromatin state distribution
#'       (\code{<mk>_chromatin_state_dist.csv})
#'     \item Sentinel: \code{.done} on successful completion
#'   }
#'
#' @details
#' Methylation profiles are always centered because 5mC/CXXC signal is
#' symmetric around CpG-dense regions. CXXC files are pulled from \code{bw_df}
#' alongside \code{mk} so both tracks appear in the same plot.
#'
#' Designed for parallel execution via \code{\link{dispatch_chromhmm_jobs}}.
#'
#' @seealso \code{\link{run_bw_profile}}, \code{\link{run_chromhmm_histone}},
#'   \code{\link{run_chromhmm_enrichment}}
#'
#' @examples
#' \dontrun{
#' bw_df <- data.frame(
#'   marker    = c("5mC", "5mC", "CXXC", "CXXC"),
#'   sample_id = c("Sample_A", "Sample_B", "Sample_A", "Sample_B"),
#'   replicate = c("rep1", "rep1", "rep1", "rep1"),
#'   batch     = c("B1", "B1", "B1", "B1"),
#'   bw_file   = c("f1.bw", "f2.bw", "f3.bw", "f4.bw")
#' )
#' loci <- GenomicRanges::GRanges("chr1:1-1000:+")
#'
#' run_chromhmm_methylation(
#'   bw_df               = bw_df,
#'   bigwig_dir          = "/path/to/bigwigs",
#'   mk                  = "5mC",
#'   loci                = loci,
#'   output_dir          = "/path/to/output",
#'   chromHmm_path       = "/path/to/chromhmm",
#'   chromHMM_annotation = "E107_15_coreMarks_hg38lift_mnemonics.bed",
#'   product             = "chromatin"
#' )
#' }
#'
#' @export
run_chromhmm_methylation <- function(bw_df, bigwig_dir, mk, loci, output_dir,
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

  # --- enrichment profile (always centered for methylation) ---
  run_bw_profile(
    allfiles      = allfiles,
    allfiles_name = allfiles_name,
    loci          = loci,
    mk            = mk,
    output_dir    = output_dir,
    mode          = "center",
    loci_label    = "CpG islands"
  )

  # --- ChromHMM state distribution (data only, no boxplot) ---
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
#' Manages parallel execution of enrichment profile and ChromHMM analysis jobs
#' using background R processes via \code{callr}. Maintains a worker pool,
#' dispatches pending jobs, and monitors completion.
#'
#' @param jobs List of job specifications, each containing:
#'   \itemize{
#'     \item \code{fn}: Function to execute (e.g.
#'       \code{run_chromhmm_histone}).
#'     \item \code{args}: Named list of arguments to pass to \code{fn}.
#'     \item \code{mk}: Character; marker name (for logging).
#'   }
#' @param n_workers Integer; maximum number of parallel workers.
#'
#' @return Invisibly returns \code{NULL} after all jobs complete.
#'
#' @details
#' \enumerate{
#'   \item Maintains a pool of up to \code{n_workers} background processes.
#'   \item Polls running processes; upon completion, checks exit status and
#'     logs results.
#'   \item Dispatches pending jobs to idle workers.
#'   \item Warns on non-zero exit; messages on success.
#'   \item Loops until all jobs are complete.
#' }
#'
#' @seealso \code{\link{run_chromhmm_enrichment}}
#'
#' @examples
#' \dontrun{
#' jobs <- list(
#'   list(
#'     fn   = run_chromhmm_histone,
#'     args = list(
#'       bw_df               = bw_df,
#'       bigwig_dir          = "/bw",
#'       mk                  = "H3K4me3",
#'       loci                = loci,
#'       output_dir          = "/out/H3K4me3",
#'       chromHmm_path       = "/path/to/chromhmm",
#'       chromHMM_annotation = "annotation.bed",
#'       product             = "chromatin"
#'     ),
#'     mk = "H3K4me3"
#'   )
#' )
#' dispatch_chromhmm_jobs(jobs, n_workers = 4)
#' }
#'
#' @export
dispatch_chromhmm_jobs <- function(jobs, n_workers) {
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
      job      <- pending[[1]]
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


#' Run ChromHMM enrichment for all markers against a set of loci
#'
#' Runs enrichment profile and ChromHMM state analysis for all markers in
#' \code{epk} against a single set of genomic \code{loci}. Call once per loci
#' set (e.g. once for protein-coding genes, once for CpG islands). Results are
#' written to \code{output_dir/<marker>/}, then loaded and stored in the EPK
#' object keyed by \code{basename(output_dir)}.
#'
#' @param epk EPK object to update.
#' @param bw_df Data frame of BigWig metadata (see \code{create_metadata_df}).
#' @param bigwig_dir Character; directory containing BigWig files.
#' @param loci \code{GRanges}; genomic features (e.g. protein-coding TSS
#'   regions, CpG islands).
#' @param output_dir Character; output directory for this loci set. The
#'   basename is used as the key in \code{epk$enrichment_results}.
#' @param chromHmm_path Character; path to the ChromHMM annotation directory.
#' @param chromHMM_annotation Character; ChromHMM annotation filename (e.g.
#'   \code{"E107_15_coreMarks_hg38lift_mnemonics.bed"}).
#' @param product Character; product type (\code{"cNUC"} or
#'   \code{"GenomePro"}).
#' @param run_mode Character; \code{"parallel"} (default) uses \code{callr}
#'   background processes; \code{"sequential"} runs one marker at a time.
#' @param n_workers Integer; parallel worker count. \code{0L} (default)
#'   auto-detects as \code{min(n_markers, detectCores() - 1)}.
#' @param markers_exclude Character vector of markers to skip.
#'   Defaults to \code{"INPUT"}.
#' @param methylation_markers Character vector of marker names routed to
#'   \code{\link{run_chromhmm_methylation}} (center-mode profile, CXXC
#'   included). All other markers use \code{\link{run_chromhmm_histone}}.
#'   Defaults to \code{c("5mC", "CXXC")}.
#' @param replicate_type Character; which replicates to use. One of
#'   \code{"replicate"} (default) or \code{"pooled"}.
#' @param plot_bw_profile_mode Character; anchoring mode for histone profiles.
#'   One of \code{"start"} (default), \code{"stretch"}, \code{"end"}, or
#'   \code{"center"}. Methylation markers always use \code{"center"}.
#'
#' @return The input \code{epk} with results appended:
#'   \itemize{
#'     \item \code{epk$enrichment_results$chromatin_states[[loci_name]]}
#'       — named list of chromatin state data frames, one per marker.
#'     \item \code{epk$enrichment_results$enrichment_profile[[loci_name]]}
#'       — named list of profile data frames, one per marker.
#'   }
#'   where \code{loci_name = basename(output_dir)}.
#'
#' @details
#' \strong{Sentinel-based skipping:} a marker is skipped if its output
#' directory already contains a \code{.done} file and all expected output
#' files. Delete \code{.done} to force re-computation.
#'
#' \strong{Routing:} markers in \code{methylation_markers} are dispatched to
#' \code{\link{run_chromhmm_methylation}}; all others go to
#' \code{\link{run_chromhmm_histone}}.
#'
#' @seealso \code{\link{run_bw_profile}}, \code{\link{run_chromhmm_histone}},
#'   \code{\link{run_chromhmm_methylation}}
#'
#' @examples
#' \dontrun{
#' epk <- run_chromhmm_enrichment(
#'   epk                 = epk, bw_df = bw_df, bigwig_dir = bigwig_dir,
#'   loci                = genes_coord_protein_coding,
#'   output_dir          = "results/Enrichment/protein_coding",
#'   chromHmm_path       = "data/chromHmm_annotations",
#'   chromHMM_annotation = "E107_15_coreMarks_hg38lift_mnemonics.bed",
#'   product             = "cNUC"
#' )
#' epk <- run_chromhmm_enrichment(
#'   epk                 = epk, bw_df = bw_df, bigwig_dir = bigwig_dir,
#'   loci                = cpg_islands_coord,
#'   output_dir          = "results/Enrichment/CpG_islands",
#'   chromHmm_path       = "data/chromHmm_annotations",
#'   chromHMM_annotation = "E107_15_coreMarks_hg38lift_mnemonics.bed",
#'   product             = "cNUC"
#' )
#' }
#'
#' @export
run_chromhmm_enrichment <- function(
  epk,
  bw_df,
  bigwig_dir,
  loci,
  output_dir,
  chromHmm_path,
  chromHMM_annotation,
  product,
  run_mode             = c("parallel", "sequential"),
  n_workers            = 0L,
  markers_exclude      = c("INPUT"),
  methylation_markers  = c("5mC", "CXXC"),
  replicate_type       = c("replicate", "pooled"),
  plot_bw_profile_mode = c("start", "stretch", "end", "center")
) {
  run_mode             <- match.arg(run_mode)
  replicate_type       <- match.arg(replicate_type)
  plot_bw_profile_mode <- match.arg(plot_bw_profile_mode)
  loci_name            <- basename(output_dir)

  markers_to_run <- setdiff(
    unique(epk$tables$stats_summary$marker),
    c(markers_exclude, NA)
  )

  if (length(markers_to_run) == 0) {
    stop("No markers to process after exclusions.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  .expected_files <- function(op, mk) {
    mode        <- if (mk %in% methylation_markers) "center" else plot_bw_profile_mode
    profile_png <- file.path(op, paste0(mk, "_profile_", mode, ".png"))
    state_csv   <- file.path(op, paste0(mk, "_chromatin_state_dist.csv"))
    c(profile_png, state_csv)
  }

  markers_needed <- Filter(function(mk) {
    op <- file.path(output_dir, mk)
    !file.exists(file.path(op, ".done")) ||
      !all(file.exists(.expected_files(op, mk)))
  }, markers_to_run)

  if (length(markers_needed) > 0) {
    nw <- if (n_workers == 0L) {
      max(1L, min(length(markers_needed), parallel::detectCores() - 1L))
    } else {
      as.integer(n_workers)
    }

    message(sprintf(
      "[chromHMM] %d marker(s) to compute for '%s': %s",
      length(markers_needed), loci_name,
      paste(markers_needed, collapse = ", ")
    ))

    common_args <- list(
      bw_df               = bw_df,
      bigwig_dir          = bigwig_dir,
      loci                = loci,
      chromHmm_path       = chromHmm_path,
      chromHMM_annotation = chromHMM_annotation,
      product             = product,
      replicate_type      = replicate_type
    )

    jobs <- lapply(markers_needed, function(mk) {
      op        <- file.path(output_dir, mk)
      dir.create(op, recursive = TRUE, showWarnings = FALSE)
      is_methyl <- mk %in% methylation_markers
      worker_fn <- if (is_methyl) run_chromhmm_methylation else run_chromhmm_histone
      extra     <- if (!is_methyl) list(plot_bw_profile_mode = plot_bw_profile_mode) else list()
      list(
        fn   = worker_fn,
        mk   = mk,
        args = c(list(mk = mk, output_dir = op), common_args, extra)
      )
    })

    if (run_mode == "parallel") {
      dispatch_chromhmm_jobs(jobs, n_workers = nw)
    } else {
      for (job in jobs) {
        message(sprintf("[chromHMM] sequential: %s", job$mk))
        callr::r(func = job$fn, args = job$args)
      }
    }
  } else {
    message(sprintf(
      "[chromHMM] all markers done for '%s', loading from cache.", loci_name
    ))
  }

  # ── Load results from disk ───────────────────────────────────────────────
  files_state <- list.files(
    output_dir,
    pattern   = ".*_chromatin_state_dist\\.csv$",
    full.names = TRUE,
    recursive  = TRUE
  )
  files_profile <- list.files(
    output_dir,
    pattern   = ".*_profile_[a-z]+_data\\.csv$",
    full.names = TRUE,
    recursive  = TRUE
  )

  .read_csvs <- function(files) {
    if (length(files) == 0) return(list())
    stats::setNames(
      lapply(files, read.csv, sep = ",", header = TRUE, stringsAsFactors = FALSE),
      basename(dirname(files))
    )
  }

  epk$enrichment_results$chromatin_states[[loci_name]]  <- .read_csvs(files_state)
  epk$enrichment_results$enrichment_profile[[loci_name]] <- .read_csvs(files_profile)

  invisible(epk)
}
