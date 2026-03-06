#' Toy metadata extracted from example BigWig filenames
#'
#' A dataset containing metadata parsed from toy BigWig file names.
#' Demonstrates the structure returned by \code{\link{create_metadata_df}}.
#'
#' @format A tibble with 6 rows and 9 columns:
#' \describe{
#'   \item{bw_file}{BigWig filename}
#'   \item{project_id}{Project identifier}
#'   \item{batch}{Batch identifier}
#'   \item{marker}{Epigenetic marker (H3K4me3, H3K9me3, H3K27me3, H3K27ac, 5mC, INPUT)}
#'   \item{rerun_id}{Rerun identifier}
#'   \item{sample_id}{Sample identifier}
#'   \item{replicate}{Replicate identifier (e.g., "rep1", "pooled")}
#'   \item{genome}{Genome version (e.g., "hg38")}
#'   \item{scaling}{Scaling method ("scaled" or "unscaled")}
#'   \item{matched}{Logical indicating if filename matched expected pattern}
#' }
#' @source Generated from toy BigWig filenames in \code{inst/extdata/toy_dataset/}
#' @examples
#' # Load toy metadata
#' data(toy_metadata)
#' head(toy_metadata)
#' 
#' # Assign to standard variable name to use with README examples
#' metadata <- toy_metadata
#' 
#' # Extract marker names
#' markers <- extract_marker_names(metadata$bw_file)
#' unique(markers)
"toy_metadata"

#' Toy QC statistics summary
#'
#' A dataset containing quality control statistics for the toy dataset samples.
#' Example output from sequencing/mapping QC pipeline.
#'
#' @format A data frame with QC metrics including:
#' \describe{
#'   \item{map_id}{BigWig filename identifier}
#'   \item{sample_id}{Sample identifier}
#'   \item{marker}{Epigenetic marker}
#'   \item{replicate}{Replicate identifier}
#'   \item{raw_demultiplexed}{Total raw reads after demultiplexing}
#'   \item{raw_mapped}{Total mapped reads}
#'   \item{mapq_mapped}{Reads with mapping quality filter}
#'   \item{dedup_mapped}{Deduplicated mapped reads}
#'   \item{final_mapped}{Final filtered mapped reads}
#'   \item{library_size}{Effective library size}
#'   \item{percent_duplication}{PCR duplication rate (percentage)}
#'   \item{frac_mapq_filtered}{Fraction of reads filtered by MAPQ}
#'   \item{insert_size}{Mean fragment insert size}
#' }
#' @source Generated from toy dataset QC pipeline
#' @examples
#' # Load toy QC statistics
#' data(toy_stats_summary)
#' head(toy_stats_summary)
#' 
#' # Assign to standard variable name to use with README examples
#' stats_summary <- toy_stats_summary
#' 
#' # Generate QC plot
#' plot_qc_stats(stats_summary, condition = "All", engine = "ggplot")
"toy_stats_summary"

#' Toy gene coordinates
#'
#' A small set of example gene coordinates for testing and examples.
#' Contains 5 protein-coding genes from chromosome 22.
#'
#' @format A GRanges object with 5 ranges and metadata columns:
#' \describe{
#'   \item{seqnames}{Chromosome (chr22)}
#'   \item{ranges}{Genomic coordinates (start-end)}
#'   \item{strand}{Strand orientation (+/-)}
#'   \item{gene_name}{Gene symbol}
#'   \item{gene_id}{Ensembl gene identifier}
#'   \item{gene_type}{Gene biotype (protein_coding)}
#' }
#' @source Subset of human genome annotations (GENCODE)
#' @examples
#' # Load toy gene coordinates
#' data(toy_genes)
#' toy_genes
#' 
#' # Assign to standard variable name to use with README examples
#' genes_coord_protein_coding <- toy_genes
#' 
#' # Use with wigglescout to process BigWig files
#' \dontrun{
#' toy_dir <- system.file("extdata", "toy_dataset", package = "EpigenicR")
#' bw_files <- list.files(toy_dir, pattern = "\\.bw$", full.names = TRUE)
#' bw_gr <- wigglescout::bw_loci(bw_files[1], loci = genes_coord_protein_coding)
#' }
"toy_genes"

#' Toy chromatin state enrichment results
#'
#' Example enrichment results showing the distribution of epigenetic marks
#' across chromatin states. Generated from toy BigWig files using 
#' wigglescout functions. Mimics the output structure from 
#' \code{wigglescout::plot_bw_loci_summary_heatmap}.
#'
#' @format A named list with 3 elements (one per marker):
#' \describe{
#'   \item{H3K4me3}{Tibble with chromatin state distribution for H3K4me3}
#'   \item{H3K27ac}{Tibble with chromatin state distribution for H3K27ac}
#'   \item{5mC}{Tibble with chromatin state distribution for DNA methylation}
#' }
#' 
#' Each element is a tibble with the following columns:
#' \describe{
#'   \item{Chromatin_State}{ChromHMM state annotation (15-state model)}
#'   \item{sample_id_rep}{Sample identifier with replicate info}
#'   \item{mean_rpgc_val}{Mean RPGC (reads per genomic content) value}
#'   \item{mean_rpgc_text}{Formatted text of mean RPGC value}
#'   \item{sample_id}{Sample identifier}
#'   \item{replicate}{Replicate/batch identifier}
#' }
#' @source Generated from toy dataset using wigglescout enrichment analysis
#' @seealso \code{\link{profile_results}} for enrichment profiles
#' @examples
#' # Load enrichment results
#' data(enrichment_results)
#' 
#' # View available markers
#' names(enrichment_results)
#' 
#' # Examine H3K4me3 enrichment across chromatin states
#' head(enrichment_results$H3K4me3)
#' 
#' # Plot chromatin state distribution
#' library(ggplot2)
#' ggplot(enrichment_results$H3K4me3, 
#'        aes(x = Chromatin_State, y = mean_rpgc_val, fill = sample_id)) +
#'   geom_bar(stat = "identity", position = "dodge") +
#'   theme_minimal() +
#'   theme(axis.text.x = element_text(angle = 45, hjust = 1))
"enrichment_results"

#' Toy enrichment profile results
#'
#' Example enrichment profiles showing epigenetic mark distribution around
#' genomic features (e.g., TSS, CpG islands). Generated from toy BigWig files 
#' using wigglescout functions. Mimics the output structure from 
#' \code{wigglescout::plot_bw_profile}.
#'
#' @format A nested named list structure:
#' \describe{
#'   \item{H3K4me3}{List with profile data for H3K4me3}
#'   \item{$protein_coding}{Profile around protein-coding gene TSS}
#'   \item{H3K27ac}{List with profile data for H3K27ac}
#'   \item{$protein_coding}{Profile around protein-coding gene TSS}
#'   \item{5mC}{List with profile data for DNA methylation}
#'   \item{$cpg_islands}{Profile around CpG islands}
#' }
#' 
#' Each profile tibble contains:
#' \describe{
#'   \item{mean_rpgc_val}{Mean RPGC signal value}
#'   \item{sderror_rpgc_val}{Standard error of RPGC}
#'   \item{median}{Median signal value}
#'   \item{index}{Distance from feature center (bp)}
#'   \item{sample}{Sample identifier with full name}
#'   \item{min_error}{Minimum error bound}
#'   \item{max_error}{Maximum error bound}
#'   \item{marker}{Epigenetic marker name}
#'   \item{sample_id}{Sample identifier}
#'   \item{sample_id_batch}{Sample identifier with batch}
#' }
#' @source Generated from toy dataset using wigglescout profile analysis
#' @seealso \code{\link{enrichment_results}} for chromatin state distributions
#' @examples
#' # Load profile results
#' data(profile_results)
#' 
#' # View available markers and feature types
#' names(profile_results)
#' names(profile_results$H3K4me3)
#' 
#' # Examine H3K4me3 profile around TSS
#' head(profile_results$H3K4me3$protein_coding)
#' 
#' # Plot enrichment profile
#' library(ggplot2)
#' ggplot(profile_results$H3K4me3$protein_coding,
#'        aes(x = index, y = mean_rpgc_val, color = sample_id)) +
#'   geom_line() +
#'   labs(x = "Distance from TSS (bp)", y = "Mean RPGC") +
#'   theme_minimal()
"profile_results"
