# Internal Setup Guide for EpigenicR

This guide is for Epigenica employees and authorized contractors setting up EpigenicR from our private GitHub repository.

## About This Package

EpigenicR is licensed under **GPL-3** (GNU General Public License version 3), ensuring:
- Any derivative works must also be open-source under GPL-3
- Full compatibility with other GPL and open-source dependencies
- Protection through access control (private GitHub repository)
- Copyright retained by Epigenica

## Installation for Epigenica Team Members

### Prerequisites

1. **GitHub Access**: Ensure you have access to Epigenica's GitHub organization
2. **Personal Access Token**: Create a GitHub PAT with `repo` scope at https://github.com/settings/tokens

### Option 1: Install from Private GitHub (Recommended)

```r
# Install from private GitHub repository
remotes::install_github("epigenica/EpigenicR", auth_token = "ghp_yourTokenHere")

# Or set up your token in .Renviron for convenience:
# GITHUB_PAT=ghp_yourTokenHere
# Then simply:
remotes::install_github("epigenica/EpigenicR")
```

### Option 2: Clone and Install Locally

```bash
# Clone the repository
git clone https://github.com/epigenica/EpigenicR.git
cd EpigenicR
```

```r
# Install from local clone
devtools::install()
```

## Setting Up Development Environment

For developers working on the package:

```r
# Load package for development
devtools::load_all("/path/to/EpigenicR")

# Run checks
devtools::check()

# Build documentation
devtools::document()

# Run tests
devtools::test()
```

## Getting Help

### Internal Resources

- **GitHub Repository**: https://github.com/epigenica/EpigenicR
- **Package Documentation**: Run `?EpigenicR` after loading
- **Example Scripts**: `inst/scripts/create_epk_example.R`
- **Support Email**: epigenica-support@epigenica.com

### Reporting Issues

For bugs or feature requests:
1. Check existing issues: https://github.com/epigenica/EpigenicR/issues
2. Create new issue with reproducible example
3. Tag with appropriate labels (bug, enhancement, etc.)
4. Mention @maintainer for urgent issues

### Contributing

1. Fork or create a branch from `main`
2. Make your changes
3. Submit a pull request
4. Wait for code review

## Data Security

### Working with Sensitive Data

- **Do not** commit actual patient/client data to version control
- Use only de-identified test datasets for examples
- Store EPK objects in secure, approved locations only
- Follow Epigenica's data handling policies

### Approved Storage Locations

- Internal server: `/epigenica/shared/analysis/`
- Personal secure storage: Check with IT for approved locations
- **NOT APPROVED**: Public cloud storage, personal devices without encryption

## Common Setup Issues

### Issue: wigglescout installation fails

```r
# Ensure correct future version first
packageurl <- "http://cran.r-project.org/src/contrib/Archive/future/future_1.34.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

# Then install wigglescout
remotes::install_github('cnluzon/wigglescout', build_vignettes = TRUE, force = TRUE)
```

### Issue: Missing Bioconductor packages

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "SummarizedExperiment", "MultiAssayExperiment"))
```

## Version History

- **v0.1.0** (March 2026): Initial internal release
  - Core functions for EPK creation
  - QC plotting
  - Annotation retrieval
  - Sample correlation analysis

## Future Development

Contact the development team if you need:
- Additional features
- Support for new genomic features
- Integration with other internal tools
- Custom analysis workflows

---

**Last Updated**: March 4, 2026  
**Maintainer**: Nimra Ali (nimra@epigenica.com)
