#!/usr/bin/env Rscript
# Test create_epk function by checking structure and dependencies

cat("=== Testing create_epk Function Structure ===\n\n")

# Read the create_epk.R file and check for syntax
source_file <- "/Users/nimra236/Dropbox/Epigenica/Projects/EpigenicR/R/create_epk.R"

cat("Reading create_epk.R...\n")
tryCatch({
  source(source_file)
  cat("✓ File loaded without syntax errors\n\n")
  
  # Check function exists
  if (exists("create_epk")) {
    cat("✓ create_epk function defined\n")
  }
  
  # Check helper functions exist
  helpers <- c(".discover_bigwig_files", ".discover_stats_summary", 
               ".normalize_annotation_to_granges", ".extract_bw_signal",
               ".extract_marker_metadata", ".get_feature_ids")
  
  for (h in helpers) {
    if (exists(h)) {
      cat("✓", h, "defined\n")
    } else {
      cat("✗", h, "NOT FOUND\n")
    }
  }
  
  cat("\n=== Function Signature ===\n")
  cat(capture.output(args(create_epk)), sep = "\n")
  
}, error = function(e) {
  cat("✗ Error loading file:\n")
  cat(conditionMessage(e), "\n")
})

cat("\n\n=== Dependencies Check ===\n")
deps <- c("wigglescout", "MultiAssayExperiment", "S4Vectors", "SummarizedExperiment", 
          "GenomicRanges", "rtracklayer")

for (dep in deps) {
  status <- ifelse(requireNamespace(dep, quietly = TRUE), "✓", "✗")
  cat(status, dep, "\n")
}
