# Example: Using the Toy Dataset
#
# This script demonstrates how to use the EpigenicR toy dataset for testing and learning.
# The toy data is loaded and assigned to standard variable names so you can use
# code examples from the README directly.

library(EpigenicR)
library(dplyr)

# ---- 1. Access toy BigWig files ----
toy_dir <- system.file("extdata", "toy_dataset", package = "EpigenicR")
bw_files <- list.files(
  file.path(toy_dir, "minute_output", "bigwig"),
  pattern = "\\.bw$",
  full.names = TRUE
)

cat("Toy dataset location:", toy_dir, "\n")
cat("Number of BigWig files:", length(bw_files), "\n\n")

# ---- 2. Load pre-computed toy data objects ----
data(toy_metadata)
data(toy_stats_summary)
data(toy_genes)

# Assign to standard variable names for easy use with README examples
metadata <- toy_metadata
stats_summary <- toy_stats_summary
genes_coord_protein_coding <- toy_genes

cat("Metadata:\n")
print(metadata)

cat("\nQC Statistics Summary:\n")
print(head(stats_summary))

cat("\nGene Coordinates (chr22 genes for testing):\n")
print(genes_coord_protein_coding)

# ---- 3. Generate QC plots ----
cat("\n---- Generating QC Plots ----\n")

# Static ggplot2 version
qc_plot_static <- plot_qc_stats(
  data = stats_summary,
  condition = "All",
  stats = c("final_mapped", "library_size", "percent_duplication"),
  marker_levels = unique(stats_summary$marker),
  engine = "ggplot",
  legend_position = "bottom",
  ncol = 3
)

print(qc_plot_static)

# Interactive plotly version
qc_plot_interactive <- plot_qc_stats(
  data = stats_summary,
  condition = "All",
  stats = c("final_mapped", "library_size"),
  engine = "plotly",
  legend_position = "bottom",
  ncol = 2
)

qc_plot_interactive

# ---- 4. Extract marker names ----
cat("\n---- Extracting Marker Names ----\n")

markers <- extract_marker_names(metadata$bw_file)
cat("Extracted markers from BigWig filenames:\n")
print(data.frame(
  sample_id = metadata$sample_id,
  marker = markers,
  bw_file = basename(metadata$bw_file)
))

cat("\nUnique markers available:", paste(unique(markers), collapse = ", "), "\n")

# ---- 5. Summary statistics by marker ----
cat("\n---- Summary Statistics by Marker ----\n")

marker_summary <- stats_summary %>%
  group_by(marker) %>%
  summarise(
    n_samples = n(),
    mean_final_reads = mean(final_mapped),
    mean_lib_size = mean(library_size),
    mean_duplication_pct = mean(percent_duplication),
    mean_insert_size = mean(insert_size)
  )

print(marker_summary)

# ---- 6. Work with specific markers ----
cat("\n---- Filtering by Specific Markers ----\n")

# Methylation samples only
methylation_samples <- metadata %>%
  filter(marker == "5mC")
cat("Methylation (5mC) samples:\n")
print(methylation_samples)

# H3K4me3 QC statistics
h3k4me3_stats <- stats_summary %>%
  filter(marker == "H3K4me3")
cat("\nH3K4me3 QC statistics:\n")
print(h3k4me3_stats)

# Chromatin marks (exclude INPUT and methylation)
chromatin_stats <- stats_summary %>%
  filter(!marker %in% c("INPUT", "5mC"))
cat("\nChromatin marks QC summary:\n")
print(chromatin_stats %>% 
        select(marker, sample_id, final_mapped, library_size, percent_duplication))

# ---- 7. Re-create metadata from BigWig filenames ----
cat("\n---- Testing create_metadata_df() Function ----\n")

# Test the metadata extraction function
metadata_test <- create_metadata_df(bw_files = basename(bw_files))
cat("Metadata extracted from filenames:\n")
print(metadata_test)

# Compare with pre-computed metadata
cat("\nMatches pre-computed metadata:", 
    all(metadata_test$marker == metadata$marker), "\n")

# ---- 8. Example: Gene information ----
cat("\n---- Gene Coordinates Information ----\n")
cat("Number of genes:", length(genes_coord_protein_coding), "\n")
cat("Chromosome:", unique(as.character(seqnames(genes_coord_protein_coding))), "\n")
cat("Gene names:", paste(head(genes_coord_protein_coding$gene_name, 5), collapse = ", "), "...\n")

# ---- Summary ----
cat("\n" ,rep("=", 60), "\n", sep = "")
cat("✓ Toy dataset exploration complete!\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nAvailable variables for further analysis:\n")
cat("  - bw_files: Paths to 6 BigWig files\n")
cat("  - metadata: Sample metadata extracted from filenames\n")
cat("  - stats_summary: QC statistics for all samples\n")
cat("  - genes_coord_protein_coding: GRanges with chr22 genes\n")
cat("  - markers: Vector of marker names\n\n")
cat("Next steps:\n")
cat("  - Process BigWig files with wigglescout::bw_loci()\n")
cat("  - Create EPK object (see: inst/scripts/create_epk_example.R)\n")
cat("  - Explore correlation analysis with compute_all_cor()\n\n")
