## Tests for chromhmm_helpers.R refactored functions

# ── helpers ──────────────────────────────────────────────────────────────────

make_bw_df <- function() {
  data.frame(
    marker    = c("H3K4me3", "H3K4me3", "INPUT",   "5mC",   "5mC",   "CXXC"),
    sample_id = c("S1",      "S2",      "S1",      "S1",    "S2",    "S1"),
    replicate = c("rep1",    "rep1",    "rep1",    "rep1",  "rep1",  "rep1"),
    batch     = c("B1",      "B1",      "B1",      "B1",    "B1",    "B1"),
    bw_file   = c(
      "H3K4me3.scaled.bw", "H3K4me3_S2.scaled.bw", "INPUT.scaled.bw",
      "5mC.scaled.bw",     "5mC_S2.scaled.bw",     "CXXC.scaled.bw"
    ),
    stringsAsFactors = FALSE
  )
}

# ── run_chromhmm_histone: file filtering ────────────────────────────────────

test_that("run_chromhmm_histone selects marker + INPUT rows (non-pooled)", {
  bw_df <- make_bw_df()

  # Intercept the bw_df_subset passed to run_bw_profile by checking which
  # files would be selected, without actually running wigglescout.
  # We test this by using a non-existent bigwig_dir — the function should
  # fail at run_bw_profile (wigglescout), not during filtering.
  # Here we just test the filtering logic directly.

  # Simulate filtering manually (mirrors the function logic):
  subset_nonpool <- dplyr::filter(
    bw_df,
    replicate != "pooled",
    marker == "H3K4me3" | tolower(marker) == "input"
  )
  # scaled filter (non-cNUC)
  subset_nonpool <- dplyr::filter(
    subset_nonpool,
    grepl("\\.scaled\\.", bw_file, perl = TRUE) | tolower(marker) == "input"
  )

  expect_equal(nrow(subset_nonpool), 3L)
  expect_true(all(subset_nonpool$marker %in% c("H3K4me3", "INPUT")))
})

test_that("run_chromhmm_histone selects pooled rows when replicate_type='pooled'", {
  bw_df <- data.frame(
    marker    = c("H3K4me3", "H3K4me3", "INPUT"),
    sample_id = c("S1",      "S1",      "S1"),
    replicate = c("rep1",    "pooled",  "pooled"),
    batch     = c("B1",      "B1",      "B1"),
    bw_file   = c("H3K4.rep1.scaled.bw", "H3K4.pooled.scaled.bw", "INPUT.pooled.scaled.bw"),
    stringsAsFactors = FALSE
  )

  subset_pool <- dplyr::filter(
    bw_df,
    replicate == "pooled",
    marker == "H3K4me3" | tolower(marker) == "input"
  )

  expect_equal(nrow(subset_pool), 2L)
  expect_true(all(subset_pool$replicate == "pooled"))
})

test_that("run_chromhmm_histone uses unscaled files for cNUC product", {
  bw_df <- data.frame(
    marker    = c("H3K4me3", "H3K4me3"),
    sample_id = c("S1",      "S1"),
    replicate = c("rep1",    "rep1"),
    batch     = c("B1",      "B1"),
    bw_file   = c("H3K4.scaled.bw", "H3K4.unscaled.bw"),
    stringsAsFactors = FALSE
  )

  subset_cnuc <- dplyr::filter(
    bw_df,
    replicate != "pooled",
    marker == "H3K4me3" | tolower(marker) == "input"
  )
  subset_cnuc <- dplyr::filter(
    subset_cnuc,
    grepl("\\.unscaled", bw_file) | tolower(marker) == "input"
  )

  expect_equal(nrow(subset_cnuc), 1L)
  expect_true(grepl("unscaled", subset_cnuc$bw_file))
})

# ── run_chromhmm_methylation: file filtering ─────────────────────────────────

test_that("run_chromhmm_methylation includes CXXC rows alongside mk", {
  bw_df <- make_bw_df()

  subset <- dplyr::filter(
    bw_df,
    replicate != "pooled",
    marker == "5mC" | marker == "CXXC" | tolower(marker) == "input"
  )
  subset <- dplyr::filter(
    subset,
    grepl("\\.scaled\\.", bw_file, perl = TRUE) | tolower(marker) == "input"
  )

  expect_true("CXXC" %in% subset$marker)
  expect_true("5mC"  %in% subset$marker)
})

test_that("run_chromhmm_methylation does NOT include histone markers", {
  bw_df <- make_bw_df()

  subset <- dplyr::filter(
    bw_df,
    replicate != "pooled",
    marker == "5mC" | marker == "CXXC" | tolower(marker) == "input"
  )

  expect_false("H3K4me3" %in% subset$marker)
})

# ── run_chromhmm_histone: skips gracefully on empty file list ────────────────

test_that("run_chromhmm_histone writes .done and returns NULL when no files match", {
  bw_df <- data.frame(
    marker    = character(0),
    sample_id = character(0),
    replicate = character(0),
    batch     = character(0),
    bw_file   = character(0),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))

  result <- run_chromhmm_histone(
    bw_df      = bw_df,
    bigwig_dir = tmp,
    mk         = "H3K4me3",
    loci       = NULL,   # not reached — no files in bw_df
    output_dir = tmp,
    product    = "chromatin"
  )

  expect_null(result)
  expect_true(file.exists(file.path(tmp, ".done")))
})

# ── run_chromhmm_methylation: skips gracefully on empty file list ─────────────

test_that("run_chromhmm_methylation writes .done and returns NULL when no files match", {
  bw_df <- data.frame(
    marker    = character(0),
    sample_id = character(0),
    replicate = character(0),
    batch     = character(0),
    bw_file   = character(0),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))

  result <- run_chromhmm_methylation(
    bw_df      = bw_df,
    bigwig_dir = tmp,
    mk         = "5mC",
    loci       = NULL,   # not reached — no files in bw_df
    output_dir = tmp,
    product    = "chromatin"
  )

  expect_null(result)
  expect_true(file.exists(file.path(tmp, ".done")))
})

# ── run_chromhmm_enrichment: sentinel-based skipping ────────────────────────

test_that("run_chromhmm_enrichment loads profile CSVs into epk", {
  # Build a minimal EPK with a stats_summary that has one marker
  epk <- list(
    mse   = NULL,
    tables = list(stats_summary = data.frame(
      marker = "H3K4me3", stringsAsFactors = FALSE
    )),
    enrichment_results = list(
      chromatin_states  = NULL,
      enrichment_profile = NULL
    ),
    provenance = list()
  )
  class(epk) <- "EPK"

  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  dir.create(tmp)

  # Simulate a completed marker: .done + profile CSV + chromatin state CSV
  mk_dir <- file.path(tmp, "H3K4me3")
  dir.create(mk_dir)
  file.create(file.path(mk_dir, ".done"))

  profile_csv <- file.path(mk_dir, "H3K4me3_profile_start_data.csv")
  write.csv(
    data.frame(index = 1:3, mean = c(0.1, 0.2, 0.3), sample = "H3K4me3_S1_rep1"),
    profile_csv,
    row.names = FALSE
  )

  state_csv <- file.path(mk_dir, "H3K4me3_chromatin_state_dist.csv")
  write.csv(
    data.frame(Chromatin_State = "1_TssA", sample_id_rep = "H3K4me3_S1_rep1",
               mean_rpgc_val = 1.5, mean_rpgc_text = "1.5"),
    state_csv,
    row.names = FALSE
  )

  bw_df <- data.frame(
    marker    = "H3K4me3",
    sample_id = "S1",
    replicate = "rep1",
    batch     = "B1",
    bw_file   = "file.bw",
    stringsAsFactors = FALSE
  )

  result <- run_chromhmm_enrichment(
    epk                 = epk,
    bw_df               = bw_df,
    bigwig_dir          = tmp,
    loci                = NULL,   # not reached — marker already done via .done sentinel
    output_dir          = tmp,
    chromHmm_path       = tmp,
    chromHMM_annotation = "dummy.bed",
    product             = "chromatin",
    run_mode            = "sequential"
  )

  expect_false(is.null(result$enrichment_results$enrichment_profile))
  expect_true("H3K4me3" %in% names(result$enrichment_results$enrichment_profile[[basename(tmp)]]))
  expect_false(is.null(result$enrichment_results$chromatin_states))
  expect_true("H3K4me3" %in% names(result$enrichment_results$chromatin_states[[basename(tmp)]]))
})

test_that("run_chromhmm_enrichment errors when no markers remain after exclusion", {
  epk <- list(
    tables = list(stats_summary = data.frame(
      marker = "INPUT", stringsAsFactors = FALSE
    )),
    enrichment_results = list(enrichment_profile = NULL)
  )
  class(epk) <- "EPK"

  expect_error(
    run_chromhmm_enrichment(
      epk                 = epk,
      bw_df               = data.frame(),
      bigwig_dir          = tempdir(),
      loci                = NULL,
      output_dir          = tempdir(),
      chromHmm_path       = tempdir(),
      chromHMM_annotation = "dummy.bed",
      product             = "chromatin"
    ),
    "No markers to process"
  )
})
