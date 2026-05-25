## Tests for chromhmm_helpers.R

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

# ── run_chromhmm_histone: file filtering ─────────────────────────────────────

test_that("run_chromhmm_histone selects marker + INPUT rows (non-pooled, scaled)", {
  bw_df <- make_bw_df()

  subset <- dplyr::filter(bw_df, replicate != "pooled",
                          marker == "H3K4me3" | tolower(marker) == "input")
  subset <- dplyr::filter(subset,
                          grepl("\\.scaled\\.", bw_file, perl = TRUE) | tolower(marker) == "input")

  expect_equal(nrow(subset), 3L)
  expect_true(all(subset$marker %in% c("H3K4me3", "INPUT")))
})

test_that("run_chromhmm_histone selects only pooled rows when replicate_type='pooled'", {
  bw_df <- data.frame(
    marker    = c("H3K4me3", "H3K4me3", "INPUT"),
    sample_id = c("S1",      "S1",      "S1"),
    replicate = c("rep1",    "pooled",  "pooled"),
    batch     = c("B1",      "B1",      "B1"),
    bw_file   = c("H3K4.rep1.scaled.bw", "H3K4.pooled.scaled.bw", "INPUT.pooled.scaled.bw"),
    stringsAsFactors = FALSE
  )

  subset <- dplyr::filter(bw_df, replicate == "pooled",
                          marker == "H3K4me3" | tolower(marker) == "input")

  expect_equal(nrow(subset), 2L)
  expect_true(all(subset$replicate == "pooled"))
})

test_that("run_chromhmm_histone selects unscaled files for cNUC product", {
  bw_df <- data.frame(
    marker    = c("H3K4me3", "H3K4me3"),
    sample_id = c("S1",      "S1"),
    replicate = c("rep1",    "rep1"),
    batch     = c("B1",      "B1"),
    bw_file   = c("H3K4.scaled.bw", "H3K4.unscaled.bw"),
    stringsAsFactors = FALSE
  )

  subset <- dplyr::filter(bw_df, replicate != "pooled",
                          marker == "H3K4me3" | tolower(marker) == "input")
  subset <- dplyr::filter(subset,
                          grepl("\\.unscaled", bw_file) | tolower(marker) == "input")

  expect_equal(nrow(subset), 1L)
  expect_true(grepl("unscaled", subset$bw_file))
})

test_that("run_chromhmm_histone does NOT include methylation markers", {
  bw_df <- make_bw_df()

  subset <- dplyr::filter(bw_df, replicate != "pooled",
                          marker == "H3K4me3" | tolower(marker) == "input")

  expect_false("5mC"  %in% subset$marker)
  expect_false("CXXC" %in% subset$marker)
})

# ── run_chromhmm_methylation: file filtering ──────────────────────────────────

test_that("run_chromhmm_methylation includes CXXC alongside mk", {
  bw_df <- make_bw_df()

  subset <- dplyr::filter(bw_df, replicate != "pooled",
                          marker == "5mC" | marker == "CXXC" | tolower(marker) == "input")
  subset <- dplyr::filter(subset,
                          grepl("\\.scaled\\.", bw_file, perl = TRUE) | tolower(marker) == "input")

  expect_true("CXXC" %in% subset$marker)
  expect_true("5mC"  %in% subset$marker)
})

test_that("run_chromhmm_methylation does NOT include histone markers", {
  bw_df <- make_bw_df()

  subset <- dplyr::filter(bw_df, replicate != "pooled",
                          marker == "5mC" | marker == "CXXC" | tolower(marker) == "input")

  expect_false("H3K4me3" %in% subset$marker)
})

# ── sentinel handling ─────────────────────────────────────────────────────────

test_that("run_chromhmm_histone writes .done and returns NULL when no files match", {
  bw_df <- data.frame(marker = character(0), sample_id = character(0),
                      replicate = character(0), batch = character(0),
                      bw_file = character(0), stringsAsFactors = FALSE)
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))

  result <- run_chromhmm_histone(
    bw_df               = bw_df,
    bigwig_dir          = tmp,
    mk                  = "H3K4me3",
    loci                = NULL,
    output_dir          = tmp,
    chromHmm_path       = tmp,
    chromHMM_annotation = "dummy.bed",
    product             = "chromatin"
  )

  expect_null(result)
  expect_true(file.exists(file.path(tmp, ".done")))
})

test_that("run_chromhmm_methylation writes .done and returns NULL when no files match", {
  bw_df <- data.frame(marker = character(0), sample_id = character(0),
                      replicate = character(0), batch = character(0),
                      bw_file = character(0), stringsAsFactors = FALSE)
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))

  result <- run_chromhmm_methylation(
    bw_df               = bw_df,
    bigwig_dir          = tmp,
    mk                  = "5mC",
    loci                = NULL,
    output_dir          = tmp,
    chromHmm_path       = tmp,
    chromHMM_annotation = "dummy.bed",
    product             = "chromatin"
  )

  expect_null(result)
  expect_true(file.exists(file.path(tmp, ".done")))
})

# ── label construction ────────────────────────────────────────────────────────

test_that("single-batch labels omit batch suffix", {
  bw_df <- data.frame(
    marker    = c("H3K4me3", "INPUT"),
    sample_id = c("S1",      "S1"),
    replicate = c("rep1",    "rep1"),
    batch     = c("B1",      "B1"),
    bw_file   = c("H3K4.scaled.bw", "INPUT.scaled.bw"),
    stringsAsFactors = FALSE
  )

  labels <- paste0(bw_df$marker, "_", bw_df$sample_id, "_", bw_df$replicate)
  expect_false(any(grepl("_B1$", labels)))
})

test_that("multi-batch labels include batch suffix", {
  bw_df <- data.frame(
    marker    = c("H3K4me3", "H3K4me3"),
    sample_id = c("S1",      "S2"),
    replicate = c("rep1",    "rep1"),
    batch     = c("B1",      "B2"),
    bw_file   = c("H3K4.scaled.bw", "H3K4_S2.scaled.bw"),
    stringsAsFactors = FALSE
  )

  labels <- paste0(bw_df$marker, "_", bw_df$sample_id,
                   "_", bw_df$replicate, "_", bw_df$batch)
  expect_true(all(grepl("_B[12]$", labels)))
})
