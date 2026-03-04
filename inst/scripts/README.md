# Example Scripts

This directory contains example scripts demonstrating how to use EpigenicR functions to create complete EPK objects.

## create_epk_example.R

This script demonstrates a complete workflow for:
- Processing BigWig files for multiple epigenomic markers
- Creating SummarizedExperiment objects for different genomic features
- Building a MultiAssayExperiment
- Bundling enrichment results
- Creating and saving an EPK object

**Note**: This is an example script that requires:
- BigWig files (`bw_files`)
- Genomic coordinates objects (genes, CpG islands, TSS regions, etc.)
- Enrichment analysis results

Adapt this script to your own data paths and analysis requirements.

## Usage

```r
# After installing EpigenicR, you can access this script:
script_path <- system.file("scripts", "create_epk_example.R", package = "EpigenicR")
file.edit(script_path)  # View/edit the script

# Copy to your working directory and adapt:
file.copy(script_path, "my_epk_creation.R")
```
