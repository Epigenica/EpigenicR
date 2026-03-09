# Using the `create_epk()` Wrapper Function

The `create_epk()` function provides two flexible input modes for creating EPK objects without manual file parsing and MultiAssayExperiment construction.

## Quick Start

### Path-Based Mode (Simplest)
If your files are organized in a pipeline output directory:

```r
library(EpigenicR)

epk <- create_epk(
  pipeline_output_path = "/path/to/pipeline/output",
  annotations = "/path/to/genes.bed"
)
```

### Explicit Mode (Full Control)
When you have individual files and want to specify exactly what goes in:

```r
epk <- create_epk(
  bw_files = c("sample1.bw", "sample2.bw", ...),
  annotations = genes_granges,
  stats_summary = qc_df
)
```

---

## Mode 1: Path-Based (Auto-Discovery)

Use this when your pipeline output follows standard directory structure:

```
project_root/
└── minute_output/
  ├── bigwig/
  │   ├── *.bw files
  └── reports/
    └── stats_summary.txt
```

### Example
```r
epk <- create_epk(
  pipeline_output_path = "/Users/research/myproject",
  annotations = "/Users/research/myproject/annotations/hg38.genes.bed"
)
```

**What it does automatically:**
- Finds all `.bw` files in `minute_output/bigwig/`
- Reads `minute_output/reports/stats_summary.txt`
- Verifies every BigWig maps to a `stats_summary$map_id` value
- Extracts marker names and sample IDs from filenames
- Creates one EPK object ready to use

---

## Mode 2: Explicit (Individual Files)

Use this for maximum control or non-standard directory layouts.

### 2a. With GRanges Annotation

```r
library(EpigenicR)

# Load from package
data(toy_genes)

# Specify files
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

# Create EPK
epk <- create_epk(
  bw_files = bw_files,
  annotations = toy_genes,      # GRanges object
  stats_summary = stats_summary
)
```

### 2b. With BED File Annotation

```r
epk <- create_epk(
  bw_files = bw_files,
  annotations = "/path/to/genes.hg38.bed",  # AutomPly read via rtracklayer
  stats_summary = stats_summary
)
```

### 2c. With Data Frame Annotation

```r
# Create annotation data frame
my_regions <- data.frame(
  chr = c("chr1", "chr1", "chr2"),
  start = c(1000, 5000, 10000),
  end = c(2000, 6000, 11000),
  feature_id = c("Gene_A", "Gene_B", "Gene_C")
)

epk <- create_epk(
  bw_files = bw_files,
  annotations = my_regions,     # Auto-converted to GRanges
  stats_summary = stats_summary
)
```

### 2d. With Multiple Annotation Sets

```r
# Create an EPK with genes AND CpG islands
epk <- create_epk(
  bw_files = bw_files,
  annotations = list(
    genes = toy_genes,
    cpg_islands = cpg_islands_granges,
    enhancers = enhancers_bed_path
  ),
  stats_summary = stats_summary
)

# Access individual experiments
experiments(epk$mse)
#> ExperimentList class object of length 3 
#> [1] "genes" "cpg_islands" "enhancers"

epk$mse[["genes"]]       # SummarizedExperiment for genes
epk$mse[["cpg_islands"]] # SummarizedExperiment for CpG islands
```

---

## Accessing Results

### Basic Info
```r
print(epk)
#> EPK object
#> ...

names(epk)
#> [1] "mse" "tables" "enrichment_results" "provenance"
```

### MultiAssayExperiment
```r
mse <- epk$mse

# Get experiments
experiments(mse)
#> ExperimentList class object of length 1 
#> [1] "primary_annotation"

# Get one experiment
se <- mse[["primary_annotation"]]
dim(se)  # rows = regions, cols = samples
#> [1] 5 6 (5 genes, 6 samples)
```

### Marker Assays
```r
se <- mse[["primary_annotation"]]

# List markers
SummarizedExperiment::assayNames(se)
#> [1] "H3K4me3" "H3K27ac" "5mC"

# Get a marker's signal matrix
h3k4me3_matrix <- SummarizedExperiment::assay(se, "H3K4me3")
dim(h3k4me3_matrix)
#> [1] 5 6
```

### QC Statistics
```r
qc_stats <- epk$tables$stats_summary
head(qc_stats)
```

### Provenance
```r
epk$provenance$created
#> [1] "2026-03-09 14:23:45 EST"

epk$provenance$session  # Full sessionInfo() at creation time
```

---

## Choose Your Mode

| Scenario | Mode | Command |
|----------|------|---------|
| Files in standard pipeline output dir | **Path** | `create_epk(pipeline_output_path = "...", annotations = "...")` |
| Individual files, simple case | **Explicit** | `create_epk(bw_files = c(...), annotations = gr, stats_summary = df)` |
| Multiple annotation sets | **Explicit** | `create_epk(bw_files = ..., annotations = list(genes = gr1, cpg = gr2), ...)` |
| Mix of BED, GRanges, data frames | **Explicit** | `create_epk(bw_files = ..., annotations = list(...), ...)` |

---

## Full Working Example

See [inst/scripts/create_epk_examples.R](../scripts/create_epk_examples.R) for runnable code with all modes.

