# Installation and Setup Guide

## Prerequisites

Before installing EpigenicR, you need to install several dependencies in the correct order.

### Step 1: Install CRAN packages

```r
install.packages(c(
  "dplyr", "tidyr", "stringr", "tibble",
  "ggplot2", "plotly", "patchwork",
  "purrr", "furrr",
  "R.utils"
))
```

### Step 2: Install Bioconductor packages

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "S4Vectors",
  "GenomicRanges",
  "SummarizedExperiment",
  "MultiAssayExperiment",
  "limma",
  "apeglm",
  "S4Arrays"
))
```

### Step 3: Install specific version of future (required for wigglescout)

```r
# Install future version 1.34.0
packageurl <- "http://cran.r-project.org/src/contrib/Archive/future/future_1.34.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
```

### Step 4: Install wigglescout from GitHub

```r
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github('cnluzon/wigglescout', build_vignettes = TRUE, force = TRUE)
```

### Step 5: Install optional packages

```r
install.packages(c("DT", "reactable", "htmltools", "ComplexHeatmap"))
```

## Installing EpigenicR

### For Development (Recommended during testing)

Use `devtools::load_all()` to load functions without installing:

```r
library(devtools)
load_all("/Users/nimra236/Dropbox/Epigenica/Projects/EpigenicR")

# Test with toy data
data(toy_genes)
epk <- create_epk(
  pipeline_output_path = system.file("extdata", "toy_dataset", package = "EpigenicR"),
  annotations = toy_genes
)
```

### From local source (Persistent installation)

```r
# Install devtools if needed
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

# Install from package directory
devtools::install_local("/Users/nimra236/Dropbox/Epigenica/Projects/EpigenicR")

# Restart R session before using
library(EpigenicR)
```

### Building and testing

After installing all dependencies:

```r
setwd("/Users/nimra236/Dropbox/Epigenica/Projects/EpigenicR")

# Generate documentation from roxygen comments
devtools::document()

# Check for errors
devtools::check()

# Run tests
devtools::test()
```

Or use devtools:

```r
devtools::document("/Users/nimra236/Dropbox/Epigenica/Projects/EpigenicR")
```

## Verifying Installation

```r
library(EpigenicR)

# Check available functions
ls("package:EpigenicR")

# View help for a function
?plot_qc_stats
```

## Troubleshooting

### Issue: wigglescout installation fails

Make sure you have the correct version of `future` (1.34.0) installed before attempting to install wigglescout.

### Issue: Bioconductor packages fail to install

Make sure you're using a compatible R version (>= 4.1.0) and have BiocManager installed.

### Issue: Documentation not generated

Ensure all package dependencies are installed before running `roxygen2::roxygenise()`.

## Development Workflow

For package development:

```r
# Load all functions without installing
devtools::load_all("/Users/nimra236/Dropbox/Epigenica/Projects/EpigenicR")

# Check package
devtools::check("/Users/nimra236/Dropbox/Epigenica/Projects/EpigenicR")

# Build package
devtools::build("/Users/nimra236/Dropbox/Epigenica/Projects/EpigenicR")
```

## Next Steps

After installation, see [README.md](README.md) for usage examples and workflows.
