# Define package root and output directories
pkg_root <- rprojroot::find_package_root_file()
output_dir <- file.path(pkg_root, "data")

# Create directories if they don't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize lists to store results
enrichment_results <- list()
profile_results <- list()

# ============================================================================
# 1. CHROMATIN STATE DISTRIBUTION (H3K4me3)
# ============================================================================

mk <- "H3K4me3"
bw_files_subset <- bw_files %>% grep(pattern = mk, value = TRUE)

bw_files_sel_name <- basename(bw_files_subset) %>%
  gsub(pattern = "\\.bw$", replacement = "") %>%
  gsub(pattern = "_pooled.*", replacement = "")

p_heatmap <- wigglescout::plot_bw_loci_summary_heatmap(
  bw_files_subset,
  "~/Epigenica/Data/chromHmm_annotation_files/E107_15_coreMarks_hg38lift_mnemonics.bed",
  labels = bw_files_sel_name,
  remove_top = 0.01
)

tmp_df <- p_heatmap@data %>%
  mutate(sample_rep = str_split(variable, "_")) %>%
  mutate(
    sample_id = gsub(sapply(sample_rep, function(x) x[5]), pattern = "\\.", replacement = "_"),
    sample_id_batch = paste0(sample_id, "_", sapply(sample_rep, function(x) x[2]))
  ) %>%
  dplyr::select(-sample_rep)

colnames(tmp_df) <- c(
  "Chromatin_State", "sample_id_rep", "mean_rpgc_val", "mean_rpgc_text", "sample_id", "replicate"
)

# Store in enrichment_results list with marker as key
enrichment_results[[mk]] <- tmp_df

# ============================================================================
# 2. ENRICHMENT PROFILE - PROTEIN CODING GENES (H3K4me3)
# ============================================================================

p_enrich <- wigglescout::plot_bw_profile(
  bw_files_subset,
  loci = genes_coord_protein_coding,
  mode = "start",
  labels = bw_files_sel_name
)

tmp_df <- p_enrich@data %>%
  mutate(sample_rep = str_split(sample, "_")) %>%
  mutate(
    marker = gsub(sapply(sample_rep, function(x) x[3]), pattern = "\\.", replacement = "_"),
    sample_id = gsub(sapply(sample_rep, function(x) x[5]), pattern = "\\.", replacement = "_"),
    sample_id_batch = paste0(sample_id, "_", sapply(sample_rep, function(x) x[2]))
  ) %>%
  dplyr::select(-sample_rep)

colnames(tmp_df) <- c(
  "mean_rpgc_val", "sderror_rpgc_val", "median", "index", "sample", "min_error", "max_error", "marker", "sample_id", "sample_id_batch"
)

# Store in profile_results list with nested structure: marker > feature_type
if (!mk %in% names(profile_results)) {
  profile_results[[mk]] <- list()
}
profile_results[[mk]][["protein_coding"]] <- tmp_df

# ============================================================================
# 3. ENRICHMENT PROFILE - CpG ISLANDS (5mC)
# ============================================================================

mk <- "5mC"
bw_files_subset <- bw_files %>% grep(pattern = mk, value = TRUE)

bw_files_sel_name <- basename(bw_files_subset) %>%
  gsub(pattern = "\\.bw$", replacement = "") %>%
  gsub(pattern = "_pooled.*", replacement = "")

# Chromatin state for 5mC (if available)
p_heatmap <- wigglescout::plot_bw_loci_summary_heatmap(
  bw_files_subset,
  "~/Epigenica/Data/chromHmm_annotation_files/E107_15_coreMarks_hg38lift_mnemonics.bed",
  labels = bw_files_sel_name,
  remove_top = 0.01
)

tmp_df <- p_heatmap@data %>%
  mutate(sample_rep = str_split(variable, "_")) %>%
  mutate(
    sample_id = gsub(sapply(sample_rep, function(x) x[5]), pattern = "\\.", replacement = "_"),
    sample_id_batch = paste0(sample_id, "_", sapply(sample_rep, function(x) x[2]))
  ) %>%
  dplyr::select(-sample_rep)

colnames(tmp_df) <- c(
  "Chromatin_State", "sample_id_rep", "mean_rpgc_val", "mean_rpgc_text", "sample_id", "replicate"
)

enrichment_results[[mk]] <- tmp_df

# Profile enrichment for CpG islands
p_enrich <- wigglescout::plot_bw_profile(
  bw_files_subset,
  loci = cpg_islands,
  mode = "center",
  labels = bw_files_sel_name
)

tmp_df <- p_enrich@data %>%
  mutate(sample_rep = str_split(sample, "_")) %>%
  mutate(
    marker = gsub(sapply(sample_rep, function(x) x[3]), pattern = "\\.", replacement = "_"),
    sample_id = gsub(sapply(sample_rep, function(x) x[5]), pattern = "\\.", replacement = "_"),
    sample_id_batch = paste0(sample_id, "_", sapply(sample_rep, function(x) x[2]))
  ) %>%
  dplyr::select(-sample_rep)

colnames(tmp_df) <- c(
  "mean_rpgc_val", "sderror_rpgc_val", "median", "index", "sample", "min_error", "max_error", "marker", "sample_id", "sample_id_batch"
)

if (!mk %in% names(profile_results)) {
  profile_results[[mk]] <- list()
}
profile_results[[mk]][["cpg_islands"]] <- tmp_df

# ============================================================================
# Save as RDA files (similar to toy_stats_summary)
# ============================================================================

usethis::use_data(enrichment_results, overwrite = TRUE)
usethis::use_data(profile_results, overwrite = TRUE)

message("Saved enrichment_results.rda and profile_results.rda to data/ directory")