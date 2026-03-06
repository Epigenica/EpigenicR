# Toy Dataset for EpigenicR Examples

This directory contains a small example dataset for demonstrating EpigenicR functionality.

## Dataset Contents

### BigWig Files (6 files)

**SAMPLE-0008** (3 markers):
1. `Proj1_A1_5mC_1_SAMPLE-0008_pooled.hg38.unscaled.bw` (3.6 MB) - Methylation
2. `Proj1_A1_H3K27ac_1_SAMPLE-0008_pooled.hg38.unscaled.bw` (1.5 MB) - H3K27 acetylation
3. `Proj1_B1_H3K4me3_1_SAMPLE-0008_pooled.hg38.unscaled.bw` (1.1 MB) - H3K4 trimethylation

**SAMPLE-0054** (3 markers):
4. `Proj1_A1_5mC_1_SAMPLE-0054_pooled.hg38.unscaled.bw` (2.0 MB) - Methylation
5. `Proj1_A1_INPUT_1_SAMPLE-0054_pooled.hg38.unscaled.bw` (2.0 MB) - INPUT control
6. `Proj1_B1_H3K4me3_1_SAMPLE-0054_pooled.hg38.unscaled.bw` (2.4 MB) - H3K4 trimethylation

**Structure**: Each sample has 3 epigenetic markers measured
**Markers included**: 5mC (2x), H3K4me3 (2x), H3K27ac (1x), INPUT (1x)  
**Genome**: hg38  
**Format**: Unscaled BigWig

### QC Statistics
- `stats_summary.txt` - Tab-delimited QC statistics table with columns:
  - map_id, library, barcode
  - raw_demultiplexed, raw_mapped, mapq_mapped
  - dedup_mapped, final_mapped, library_size
  - percent_duplication, frac_mapq_filtered, insert_size

## Accessing Toy Data

After package installation, users can access these files and pre-computed R data objects:

### Quick Start (Recommended)

Load the toy data and assign to standard variable names for immediate use:

```r
library(EpigenicR)

# Get toy BigWig file paths
toy_dir <- system.file("extdata", "toy_dataset", package = "EpigenicR")
bw_files <- list.files(toy_dir, pattern = "\\.bw$", full.names = TRUE)

# Load pre-computed R data objects
data(toy_metadata)
data(toy_stats_summary)
data(toy_genes)

# Assign to standard variable names (makes README examples work directly!)
metadata <- toy_metadata
stats_summary <- toy_stats_summary
genes_coord_protein_coding <- toy_genes

# Now use any example from the main README
qc_plot <- plot_qc_stats(
  data = stats_summary,
  condition = "All",
  engine = "ggplot"
)
print(qc_plot)
```

### Accessing Individual Files

```r
# List available BigWig files
toy_dir <- system.file("extdata", "toy_dataset", package = "EpigenicR")
list.files(toy_dir, pattern = "\\.bw$")

# Get specific file path
bw_file <- system.file("extdata", "toy_dataset", 
                       "Proj1_A1_5mC_1_SAMPLE-0008_pooled.hg38.unscaled.bw", 
                       package = "EpigenicR")

# Load stats summary directly from file
stats_file <- system.file("extdata", "toy_dataset", "stats_summary.txt", 
                          package = "EpigenicR")
stats <- read.table(stats_file, header = TRUE, sep = "\t")
```

### Why Assign to Standard Names?

The toy data objects are prefixed with `toy_` for clarity in the package namespace, but this makes it harder to use code examples directly. By assigning:

```r
metadata <- toy_metadata
stats_summary <- toy_stats_summary  
genes_coord_protein_coding <- toy_genes
```

You can copy-paste any code example from the README and it will work immediately with the toy dataset!

## File Size Considerations

Keep files small for package distribution:
- BigWig files: Ideally < 5 MB each (subsample if needed)
- Consider only including a small genomic region (e.g., chr22 or specific genes)
- Total package size should be < 50 MB for CRAN-like repositories

## Usage in Examples

These files are used in:
- Package documentation examples
- Vignettes
- Test cases
- README demonstrations
