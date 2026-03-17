#' Run ChromHMM enrichment analysis for histone marks
#'
#' Computes enrichment profiles and chromatin state distributions for histone
#' modification markers (e.g., H3K4me3, H3K27ac) at genomic features using
#' ChromHMM annotations. Generates visualization and summary tables.
#'
#' @param bw_df Data frame containing BigWig file metadata with columns:
#'   \code{marker}, \code{sample_id}, \code{replicate}, \code{batch},
#'   \code{scaling}, \code{bw_file}.
#' @param bigwig_dir Character; path to directory containing BigWig files.
#' @param mk Character; marker name (e.g., "H3K4me3", "H3K27ac").
#' @param loci A \code{GRanges} object specifying genomic features for analysis.
#' @param output_dir Character; path where output plots and tables are saved.
#' @param chromHmm_path Character; path to ChromHMM annotation directory.
#' @param chromHMM_annotation Character; name of ChromHMM annotation file
#'   (e.g., "E107_15_coreMarks_hg38lift_mnemonics.bed").
#' @param product Character; product type ("cNUC" uses unscaled bigWigs;
#'   others use scaled).
#'
#' @return Invisibly returns \code{NULL}. Writes:
#'   \itemize{
#'     \item PNG: enrichment profile plot (\code{<marker>_profile_start.png})
#'     \item CSV: profile data (\code{<marker>_profile_start_data.csv})
#'     \item PNG: chromatin state distribution (\code{<marker>_chromatin_state_dist.png})
#'     \item CSV: chromatin state summary (\code{<marker>_chromatin_state_dist.csv})
#'     \item Marker file: \code{.done} upon completion
#'   }
#'
#' @details
#' Generates two outputs:
#' \enumerate{
#'   \item \strong{Enrichment profile}: Signal at TSS for pooled replicates.
#'   \item \strong{Chromatin state distribution}: Mean RPGC per state across replicates.
#' }
#'
#' Designed for parallel execution via \code{dispatch_chromhmm_jobs()}.
#'
#' @examples
#' \dontrun{
#' # Setup sample data
#' bw_df <- data.frame(
#'   marker = c("H3K4me3", "H3K4me3", "INPUT"),
#'   sample_id = c("Sample_A", "Sample_A", "Sample_A"),
#'   replicate = c("rep1", "pooled", "pooled"),
#'   batch = c("B1", "B1", "B1"),
#'   scaling = c("scaled", "scaled", "scaled"),
#'   bw_file = c("file1.bw", "file2.bw", "file3.bw")
#' )
#' loci <- GenomicRanges::GRanges("chr1:1-1000:+")
#'
#' run_chromhmm_histone_enrichment(
#'   bw_df = bw_df,
#'   bigwig_dir = "/path/to/bigwigs",
#'   mk = "H3K4me3",
#'   loci = loci,
#'   output_dir = "/path/to/output",
#'   chromHmm_path = "/path/to/chromhmm",
#'   chromHMM_annotation = "E107_15_coreMarks_hg38lift_mnemonics.bed",
#'   product = "chromatin"
#' )
#' }
#'
#' @export
run_chromhmm_histone_enrichment <- function(bw_df, bigwig_dir, mk, loci,
                                 output_dir, chromHmm_path, chromHMM_annotation,
                                 product) {
  library(GenomicRanges)
  library(wigglescout)
  library(ggplot2)
  library(dplyr)
  library(stringr)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # --- enrichment profile ---
  bw_df_subset <- bw_df %>%
    dplyr::filter(replicate == "pooled") %>%
    dplyr::filter(marker %in% c(mk, "INPUT"))

  allfiles <- file.path(bigwig_dir, bw_df_subset$bw_file)
  if (length(allfiles) == 0) {
    message(sprintf("[chromHMM:%s] no pooled bigWig files found in bw_df - skipping enrichment profile.", mk))
  } else {
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

    p_enrich <- wigglescout::plot_bw_profile(
      allfiles,
      loci = loci,
      mode = "start",
      labels = allfiles_name
    )

    p_enrich <- p_enrich +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::labs(x = paste0(mk, " Protein coding (", length(loci), ") loci")) +
      ggplot2::guides(color = ggplot2::guide_legend(ncol = 1, byrow = TRUE)) +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_blank()
      )

    ggplot2::ggsave(
      file.path(output_dir, paste0(mk, "_profile_start.png")),
      plot = p_enrich,
      width = 12,
      height = 8
    )

    write.table(
      p_enrich$data,
      file.path(output_dir, paste0(mk, "_profile_start_data.csv")),
      row.names = FALSE,
      quote = FALSE,
      col.names = TRUE,
      sep = ","
    )
  }

  # --- chromHMM boxplot ---
  scaling_type <- if (product != "cNUC") "scaled" else "unscaled"

  bw_df_subset <- bw_df %>%
    dplyr::filter(
      marker %in% c(mk, "INPUT"),
      scaling == scaling_type,
      replicate != "pooled"
    )

  bw_files_sel_path <- file.path(bigwig_dir, bw_df_subset$bw_file)
  bw_files_sel_name <- paste0(bw_df_subset$sample_id, "_", bw_df_subset$replicate)

  if (length(bw_files_sel_path) == 0) {
    message(sprintf("[chromHMM:%s] no scaled non-pooled bigWig files found in bw_df - skipping chromHMM boxplot.", mk))
  } else {
    p_heatmap <- wigglescout::plot_bw_loci_summary_heatmap(
      bw_files_sel_path,
      file.path(chromHmm_path, chromHMM_annotation),
      labels = bw_files_sel_name,
      remove_top = 0.01
    )

    tmp_df <- p_heatmap@data %>%
      dplyr::mutate(sample_rep = stringr::str_split(variable, "_")) %>%
      dplyr::mutate(
        sample_id = gsub(sapply(sample_rep, `[`, 1), pattern = "\\.", replacement = "_"),
        replicate = sapply(sample_rep, `[`, 2)
      ) %>%
      dplyr::select(-sample_rep)

    colnames(tmp_df) <- c(
      "Chromatin_State", "sample_id_rep", "mean_rpgc_val",
      "mean_rpgc_text", "sample_id", "replicate"
    )

    p_box <- ggplot2::ggplot(tmp_df, ggplot2::aes(x = Chromatin_State, y = mean_rpgc_val)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(
        title = paste0(mk, " coverage per region (mean RPGC)"),
        x = "Features",
        y = "RPGC"
      ) +
      ggplot2::facet_wrap(~sample_id, scales = "free_x") +
      ggplot2::coord_flip()

    ggplot2::ggsave(
      file.path(output_dir, paste0(mk, "_chromatin_state_dist.png")),
      plot = p_box,
      width = 12,
      height = 8
    )

    write.table(
      tmp_df,
      file.path(output_dir, paste0(mk, "_chromatin_state_dist.csv")),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
  }

}

#' Run ChromHMM enrichment analysis for methylation marks
#'
#' Computes enrichment profiles and chromatin state distributions for methylation
#' markers (5mC and CXXC) at genomic features using ChromHMM annotations.
#' Generates visualization and summary tables.
#'
#' @param bw_df Data frame containing BigWig file metadata with columns:
#'   \code{marker}, \code{sample_id}, \code{replicate}, \code{batch},
#'   \code{scaling}, \code{bw_file}.
#' @param bigwig_dir Character; path to directory containing BigWig files.
#' @param mk Character; marker name (typically "5mC" or "CXXC").
#' @param loci A \code{GRanges} object specifying genomic features for analysis.
#' @param output_dir Character; path where output plots and tables are saved.
#' @param chromHmm_path Character; path to ChromHMM annotation directory.
#' @param chromHMM_annotation Character; name of ChromHMM annotation file
#'   (e.g., "E107_15_coreMarks_hg38lift_mnemonics.bed").
#' @param product Character; product type ("cNUC" uses unscaled bigWigs;
#'   others use scaled).
#'
#' @return Invisibly returns \code{NULL}. Writes:
#'   \itemize{
#'     \item PNG: enrichment profile plot (\code{methylation_profile_center.png})
#'     \item CSV: profile data (\code{methylation_profile_center_data.csv})
#'     \item PNG: chromatin state distribution (\code{5mC_chromatin_state_dist.png})
#'     \item CSV: chromatin state summary (\code{5mC_chromatin_state_dist.csv})
#'     \item Marker file: \code{.done} upon completion
#'   }
#'
#' @details
#' Generates two outputs:
#' \enumerate{
#'   \item \strong{Enrichment profile}: Signal at central region for pooled replicates of 5mC and CXXC.
#'   \item \strong{Chromatin state distribution}: Mean RPGC per state across 5mC replicates.
#' }
#'
#' Designed for parallel execution via \code{dispatch_chromhmm_jobs()}.
#'
#' @examples
#' \dontrun{
#' # Setup sample data
#' bw_df <- data.frame(
#'   marker = c("5mC", "5mC", "CXXC", "CXXC"),
#'   sample_id = c("Sample_A", "Sample_A", "Sample_A", "Sample_A"),
#'   replicate = c("rep1", "pooled", "rep1", "pooled"),
#'   batch = c("B1", "B1", "B1", "B1"),
#'   scaling = c("scaled", "scaled", "scaled", "scaled"),
#'   bw_file = c("file1.bw", "file2.bw", "file3.bw", "file4.bw")
#' )
#' loci <- GenomicRanges::GRanges("chr1:1-1000:+")
#'
#' run_chromhmm_methylation_enrichment(
#'   bw_df = bw_df,
#'   bigwig_dir = "/path/to/bigwigs",
#'   mk = "5mC",
#'   loci = loci,
#'   output_dir = "/path/to/output",
#'   chromHmm_path = "/path/to/chromhmm",
#'   chromHMM_annotation = "E107_15_coreMarks_hg38lift_mnemonics.bed",
#'   product = "chromatin"
#' )
#' }
#'
#' @export
run_chromhmm_methylation_enrichment <- function(bw_df, bigwig_dir, mk, loci,
                                    output_dir, chromHmm_path, chromHMM_annotation,
                                    product) {
  library(GenomicRanges)
  library(wigglescout)
  library(ggplot2)
  library(dplyr)
  library(stringr)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # --- enrichment profile ---
  bw_df_subset <- bw_df %>%
    dplyr::filter(replicate == "pooled") %>%
    dplyr::filter(marker %in% c("5mC", "CXXC"))

  allfiles <- file.path(bigwig_dir, bw_df_subset$bw_file)

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

  p_enrich <- wigglescout::plot_bw_profile(
    allfiles,
    loci = loci,
    mode = "center",
    labels = allfiles_name
  )

  p_enrich <- p_enrich +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::labs(x = paste0(mk, " Protein coding (", length(loci), ") loci")) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 1, byrow = TRUE)) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(output_dir, "methylation_profile_center.png"),
    plot = p_enrich,
    width = 12,
    height = 8
  )

  write.table(
    p_enrich$data,
    file.path(output_dir, "methylation_profile_center_data.csv"),
    row.names = FALSE,
    quote = FALSE,
    col.names = TRUE,
    sep = ","
  )

  # --- chromHMM boxplot ---
  scaling_type <- if (product != "cNUC") "scaled" else "unscaled"
  bw_df_subset <- dplyr::filter(bw_df, marker %in% "5mC", scaling == scaling_type)

  if (length(unique(bw_df_subset$batch)) > 1) {
    bw_files_sel_name <- paste0(
      bw_df_subset$sample_id, "_", bw_df_subset$replicate,
      "_", bw_df_subset$batch, "_", bw_df_subset$marker
    )
  } else {
    bw_files_sel_name <- paste0(
      bw_df_subset$sample_id, "_", bw_df_subset$replicate,
      "_", bw_df_subset$marker
    )
  }

  bw_files_sel_path <- file.path(bigwig_dir, bw_df_subset$bw_file)

  p_heatmap <- wigglescout::plot_bw_loci_summary_heatmap(
    bw_files_sel_path,
    file.path(chromHmm_path, chromHMM_annotation),
    labels = bw_files_sel_name,
    remove_top = 0.01
  )

  tmp_df <- p_heatmap@data %>%
    dplyr::mutate(sample_rep = stringr::str_split(variable, "_")) %>%
    dplyr::mutate(
      sample_id = gsub(sapply(sample_rep, `[`, 1), pattern = "\\.", replacement = "_"),
      replicate = sapply(sample_rep, `[`, 2)
    ) %>%
    dplyr::select(-sample_rep)

  colnames(tmp_df) <- c(
    "Chromatin_State", "sample_id_rep", "mean_rpgc_val",
    "mean_rpgc_text", "sample_id", "replicate"
  )

  p_box <- ggplot2::ggplot(tmp_df, ggplot2::aes(x = Chromatin_State, y = mean_rpgc_val)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(
      title = paste0(mk, " coverage per region (mean RPGC)"),
      x = "Features",
      y = "RPGC"
    ) +
    ggplot2::facet_wrap(~sample_id, scales = "free_x") +
    ggplot2::coord_flip()

  ggplot2::ggsave(
    file.path(output_dir, paste0(mk, "_chromatin_state_dist.png")),
    plot = p_box,
    width = 12,
    height = 8
  )

  write.table(
    tmp_df,
    file.path(output_dir, paste0(mk, "_chromatin_state_dist.csv")),
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  file.create(file.path(output_dir, ".done"))
  invisible(NULL)
}

#' Dispatch and orchestrate parallel ChromHMM jobs
#'
#' Manages parallel execution of ChromHMM enrichment analysis jobs using
#' background R processes via \code{callr}. Maintains a worker pool, dispatches
#' pending jobs, and monitors completion.
#'
#' @param jobs List of job specifications, each containing:
#'   \itemize{
#'     \item \code{fn}: Function to execute (e.g., \code{run_chromhmm_histone_enrichment}).
#'     \item \code{args}: Named list of arguments to pass to \code{fn}.
#'     \item \code{mk}: Character; marker name (for logging).
#'   }
#' @param n_workers Integer; maximum number of parallel workers to maintain.
#'
#' @return Invisibly returns \code{NULL} after all jobs complete.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Maintains a pool of up to \code{n_workers} background processes.
#'   \item Polls running processes; upon completion, checks exit status and logs results.
#'   \item Dispatches pending jobs to idle workers.
#'   \item Warns on non-zero exit status; messages on success.
#'   \item Loops until all jobs are complete.
#' }
#'
#' Designed for large-scale, batch ChromHMM analysis where multiple markers
#' can be processed simultaneously.
#'
#' @examples
#' \dontrun{
#' # Setup job list
#' jobs <- list(
#'   list(
#'     fn = run_chromhmm_histone_enrichment,
#'     args = list(
#'       bw_df = bw_df,
#'       bigwig_dir = "/path/to/bw",
#'       mk = "H3K4me3",
#'       loci = loci,
#'       output_dir = "/path/to/out",
#'       chromHmm_path = "/path/to/chromhmm",
#'       chromHMM_annotation = "annotation.bed",
#'       product = "chromatin"
#'     ),
#'     mk = "H3K4me3"
#'   ),
#'   list(
#'     fn = run_chromhmm_histone_enrichment,
#'     args = list(...),
#'     mk = "H3K27ac"
#'   )
#' )
#'
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
      j <- running[[i]]
      err <- j$proc$read_all_error()

      if (j$proc$get_exit_status() != 0) {
        warning(
          sprintf(
            "[chromHMM] worker '%s' FAILED (exit %d):\n%s",
            j$mk,
            j$proc$get_exit_status(),
            err
          )
        )
      } else {
        if (nzchar(trimws(err))) {
          message(sprintf("[chromHMM] worker '%s' stderr:\n%s", j$mk, err))
        }
        message(sprintf("[chromHMM] completed: %s", j$mk))
      }
    }

    running <- running[!done]

    while (length(running) < n_workers && length(pending) > 0) {
      job <- pending[[1]]
      pending[[1]] <- NULL

      proc <- callr::r_bg(func = job$fn, args = job$args)
      running <- c(running, list(list(proc = proc, mk = job$mk)))

      message(sprintf(
        "[chromHMM] dispatched: %s (%d/%d running)",
        job$mk,
        length(running),
        n_workers
      ))
    }

    if (length(running) == 0 && length(pending) == 0) {
      break
    }

    Sys.sleep(1)
  }

  invisible(NULL)
}
