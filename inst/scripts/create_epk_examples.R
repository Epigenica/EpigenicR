# Examples: Using create_epk() Wrapper Function
# 
# This script demonstrates both input modes for the create_epk() wrapper:
# 1. Path-based mode (auto-discovery)
# 2. Explicit mode (individual files)

# ==============================================================================
# MODE 1: PATH-BASED (Auto-discovery from pipeline output)
# ==============================================================================

library(EpigenicR)

# Example: Running on real pipeline output structure
# Assume your pipeline output is organized like:
#   /path/to/project/
#   └── minute_output/
#       ├── bigwig/
#       │   ├── Proj1_A1_H3K4me3_1_SAMPLE-0008_rep1.hg38.unscaled.bw
#       │   ├── Proj1_A1_H3K27ac_1_SAMPLE-0008_rep1.hg38.unscaled.bw
#       │   └── ... (more .bw files)
#       └── reports/
#           └── stats_summary.txt

# One-liner: Create EPK from pipeline output with annotation file
epk_from_pipeline <- create_epk(
  pipeline_output_path = "/path/to/project",
  annotations = "/path/to/regions.bed"
)

# Access results
print(epk_from_pipeline)
MultiAssayExperiment::experiments(epk_from_pipeline$mse)


# ==============================================================================
# MODE 2: EXPLICIT (Individual files)
# ==============================================================================

# Example 2a: Using toy dataset from package
data(toy_genes)

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

# Create EPK with explicit files and GRanges annotation
epk_explicit <- create_epk(
  bw_files = bw_files,
  annotations = toy_genes,
  stats_summary = stats_summary
)

print(epk_explicit)


# ==============================================================================
# MODE 2b: EXPLICIT with BED file annotation
# ==============================================================================

epk_explicit_bed <- create_epk(
  bw_files = bw_files,
  annotations = "/path/to/annotation.bed",
  stats_summary = stats_summary
)


# ==============================================================================
# MODE 2c: EXPLICIT with data.frame annotation
# ==============================================================================

# Create a simple annotation data frame
annotation_df <- data.frame(
  chr = c("chr1", "chr1", "chr2"),
  start = c(1000, 5000, 10000),
  end = c(2000, 6000, 11000),
  feature_id = c("region_A", "region_B", "region_C")
)

epk_explicit_df <- create_epk(
  bw_files = bw_files,
  annotations = annotation_df,
  stats_summary = stats_summary
)


# ==============================================================================
# MODE 2d: EXPLICIT with multiple annotation sets
# ==============================================================================

# Load additional annotations
data(toy_genes)

# Create BED file paths or load from disk
genes_annotation <- toy_genes
# cpg_annotation <- rtracklayer::import("path/to/cpg_islands.bed")

# Create EPK with multiple experiments
epk_multi <- create_epk(
  bw_files = bw_files,
  annotations = list(
    genes = toy_genes
    # cpg_islands = cpg_annotation
  ),
  stats_summary = stats_summary
)

# Access individual experiments
experiments(epk_multi$mse)
epk_multi$mse[["genes"]]  # SummarizedExperiment for genes


# ==============================================================================
# POST-CREATION: Access EPK components
# ==============================================================================

# View EPK structure
print(epk_explicit)

# Access MultiAssayExperiment
mse <- epk_explicit$mse

# Access QC statistics
qc_stats <- epk_explicit$tables$stats_summary
head(qc_stats)

# Access individual experiments
se_genes <- MultiAssayExperiment::experiments(mse)[[1]]
dim(se_genes)  # rows = regions, cols = samples

# Get assay matrices (one per marker)
markernames <- SummarizedExperiment::assayNames(se_genes)
print(markernames)

# Extract H3K4me3 signal matrix
h3k4me3_matrix <- SummarizedExperiment::assay(se_genes, "H3K4me3")
print(dim(h3k4me3_matrix))

# Compute correlations between samples
cor_h3k4me3 <- cor(h3k4me3_matrix, method = "pearson")
print(cor_h3k4me3)
