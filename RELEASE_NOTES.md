# EpigenicR v0.1.0 - Initial Release 🎉

First official release of EpigenicR for epigenomic data analysis from the EpiFinder platform!

## 🚀 Key Features

### Data Processing
- **Metadata extraction** from BigWig filenames with automatic parsing
- **QC visualization** with interactive and static plotting options
- **ChromHMM integration** for chromatin state annotations
- **Genomic annotation** management (GTF/BED files)

### Data Structures
- **EPK objects** containing MultiAssayExperiment with enrichment results
- Support for multiple markers: H3K4me3, H3K27ac, H3K9me3, H3K27me3, 5mC, INPUT
- Integration with wigglescout for BigWig processing

### Analysis Tools
- Sample-sample **correlation analysis** across markers
- Chromatin state **enrichment profiles**
- Signal distribution around genomic features (TSS, CpG islands, enhancers)

## 📦 Complete Example Dataset

Ready-to-use toy dataset with:
- 8 real BigWig files (H3K4me3, H3K27ac, 5mC, INPUT from 2 samples)
- QC statistics and metadata
- Pre-computed enrichment results
- Example EPK object
- All loadable with `data()`

## 📚 Documentation

- Comprehensive README with workflows
- Complete function documentation
- Getting started vignette
- Example analysis scripts

## 📥 Installation

```r
# Install dependencies
packageurl <- "http://cran.r-project.org/src/contrib/Archive/future/future_1.34.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
remotes::install_github('cnluzon/wigglescout', build_vignettes = TRUE, force = TRUE)

# Install EpigenicR
remotes::install_github("epigenica/EpigenicR@v0.1.0")
```

## 🎯 Quick Start

```r
library(EpigenicR)

# Load toy data
data(toy_metadata)
data(toy_stats_summary)
data(enrichment_results)

# Visualize QC
plot_qc_stats(toy_stats_summary, condition = "All", engine = "ggplot")

# View enrichment results
names(enrichment_results)  # H3K4me3, H3K27ac, 5mC
```

## 📋 What's Included

### Core Functions
- `create_metadata_df()` - Parse BigWig filenames
- `plot_qc_stats()` - Generate QC plots
- `ensure_gtf_and_beds()` - Download genomic annotations
- `download_chromhmm_annotations()` - Get chromatin states
- `compute_sample_cor()` - Calculate correlations
- `extract_marker_names()` - Extract marker info

### Example Data
- `toy_metadata` - Metadata from 6 samples
- `toy_stats_summary` - QC statistics for 8 samples
- `toy_genes` - Example gene coordinates (chr22)
- `enrichment_results` - Chromatin state distributions
- `profile_results` - Enrichment profiles
- `project.epk.rds` - Complete EPK object example

### Scripts
- `create_epk_example.R` - Build EPK objects
- `wigglescout_script.R` - Enrichment analysis workflow
- `generate_example_plots.R` - Plotting examples
- `toy_dataset_examples.R` - Work with toy data

## 🔧 System Requirements

- R >= 4.1.0
- future v1.34.0 (specific version required)
- wigglescout (GitHub version)
- Standard Bioconductor packages

## 📝 Notes

- Designed for hg38 genome build
- Optimized for EpiFinder platform data
- Private repository for internal Epigenica use

---

**Full Changelog**: See [NEWS.md](NEWS.md) for detailed changes
