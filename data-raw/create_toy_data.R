# Process Toy Dataset
# This script creates R data objects from the toy dataset files
# Run this after uploading the BigWig files and stats_summary.txt

# Load required packages
library(GenomicRanges)
library(dplyr)
library(tibble)

# Source functions directly (for development)
source("R/functions.R")

# ---- Load toy data paths ----
toy_dir <- "inst/extdata/toy_dataset"
toy_minute_output_dir <- file.path(toy_dir, "minute_output")
toy_bigwig_dir <- file.path(toy_minute_output_dir, "bigwig")
toy_stats_file <- file.path(toy_minute_output_dir, "reports", "stats_summary.txt")

# ---- Create toy_bw_files object ----
toy_bw_files <- list.files(toy_bigwig_dir, pattern = "\\.bw$", full.names = TRUE)
names(toy_bw_files) <- basename(toy_bw_files)

message("Found ", length(toy_bw_files), " BigWig files:")
print(basename(toy_bw_files))

# ---- Create toy_metadata ----
toy_metadata <- create_metadata_df(bw_files = basename(toy_bw_files))
message("\nGenerated metadata:")
print(toy_metadata)

# ---- Load toy stats_summary ----
toy_stats_summary <- read.table(
  toy_stats_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Add marker column for easier plotting
toy_stats_summary <- toy_stats_summary %>%
  mutate(
    marker = case_when(
      grepl("5mC", map_id) ~ "5mC",
      grepl("H3K4me3", map_id) ~ "H3K4me3",
      grepl("H3K27ac", map_id) ~ "H3K27ac",
      grepl("INPUT", map_id) ~ "INPUT",
      TRUE ~ "Unknown"
    ),
    sample_id = gsub(".*_(SAMPLE-[0-9]+).*", "\\1", map_id),
    replicate = sub(".*_(rep[0-9]+)\\..*", "\\1", map_id),
    condition = "toy_dataset"
  )

message("\nProcessed stats_summary with ", nrow(toy_stats_summary), " samples")

# ---- Create small example genomic ranges ----
# Example: chr22 genes (small chromosome, good for examples)
toy_genes <- GRanges(
  seqnames = "chr22",
  ranges = IRanges(
    start = c(10500000, 10600000, 10700000, 10800000, 10900000),
    end =   c(10510000, 10610000, 10710000, 10810000, 10910000)
  ),
  strand = c("+", "-", "+", "+", "-"),
  gene_name = c("TOY1", "TOY2", "TOY3", "TOY4", "TOY5"),
  gene_id = c("ENSG00000000001", "ENSG00000000002", "ENSG00000000003", 
              "ENSG00000000004", "ENSG00000000005"),
  gene_type = rep("protein_coding", 5)
)

# ---- Save as R data objects ----
# These will be available to users as data(toy_metadata)
usethis::use_data(toy_metadata, overwrite = TRUE)
usethis::use_data(toy_stats_summary, overwrite = TRUE)
usethis::use_data(toy_genes, overwrite = TRUE)

message("\n✓ Toy data objects created successfully!")
message("Objects saved to data/ directory:")
message("  - toy_metadata (", nrow(toy_metadata), " files)")
message("  - toy_stats_summary (", nrow(toy_stats_summary), " samples)")
message("  - toy_genes (", length(toy_genes), " genes)")
message("\nUsers can load with: data(toy_metadata)")

