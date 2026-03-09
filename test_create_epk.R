#!/usr/bin/env Rscript
# Quick test of create_epk wrapper function

# Load the package
library(devtools)
load_all()

# Check that create_epk is exported
if (!exists("create_epk")) {
  stop("create_epk function not found!")
}

cat("✓ create_epk function successfully loaded\n")

# Test with toy data
cat("\n\n=== Testing create_epk with toy dataset ===\n")

data(toy_metadata, toy_genes)

toy_dir <- system.file("extdata", "toy_dataset", package = "EpigenicR")
bw_files <- list.files(
  file.path(toy_dir, "minute_output", "bigwig"),
  pattern = "\\.bw$",
  full.names = TRUE
)
stats_summary <- read.table(
  file.path(toy_dir, "minute_output", "reports", "stats_summary.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

cat("BigWig files found:", length(bw_files), "\n")
cat("Sample stats_summary:", nrow(stats_summary), "rows\n")
cat("Annotations (toy_genes):", length(toy_genes), "regions\n")

cat("\nAttempting to create EPK object...\n")

tryCatch({
  epk <- create_epk(
    bw_files = bw_files,
    annotations = toy_genes,
    stats_summary = stats_summary
  )
  
  cat("\n✓ EPK object created successfully!\n")
  cat("EPK class:", class(epk), "\n")
  cat("EPK components:", names(epk), "\n")
  cat("MultiAssayExperiment experiments:", names(MultiAssayExperiment::experiments(epk$mse)), "\n")
  
}, error = function(e) {
  cat("\n✗ Error creating EPK:\n")
  cat(conditionMessage(e), "\n")
  traceback()
})
