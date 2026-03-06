# Create Toy Enrichment and Profile Results
# This script creates example enrichment_results and profile_results data objects
# that mimic the structure of wigglescout output for the toy dataset

library(dplyr)
library(tibble)

# ---- Create toy enrichment_results ----
# Structure: list with marker names as keys, containing chromatin state distribution data

# Example chromatin states from ChromHMM (15-state model)
chromatin_states <- c(
  "1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD", "5_Tx",
  "6_TxWk", "7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2",
  "11_EnhWk", "12_ZNF/Rpts", "13_Het", "14_TssBiv", "15_EnhBiv"
)

# Create toy data for H3K4me3 (promoter marker)
h3k4me3_enrichment <- tibble(
  Chromatin_State = rep(chromatin_states, 2),
  sample_id_rep = rep(c("Proj1_A1_H3K4me3_1_SAMPLE_0008", 
                        "Proj1_A1_H3K4me3_1_SAMPLE_0054"), each = 15),
  mean_rpgc_val = c(
    # Sample 0008 - high at promoters, low elsewhere
    2.5, 1.8, 1.2, 1.0, 0.8, 0.6, 0.4, 0.3, 0.3, 0.2, 0.2, 0.1, 0.1, 0.5, 0.3,
    # Sample 0054 - similar pattern
    2.3, 1.7, 1.1, 0.9, 0.7, 0.5, 0.35, 0.25, 0.28, 0.18, 0.18, 0.12, 0.09, 0.48, 0.28
  ),
  mean_rpgc_text = sprintf("%.2f", c(
    2.5, 1.8, 1.2, 1.0, 0.8, 0.6, 0.4, 0.3, 0.3, 0.2, 0.2, 0.1, 0.1, 0.5, 0.3,
    2.3, 1.7, 1.1, 0.9, 0.7, 0.5, 0.35, 0.25, 0.28, 0.18, 0.18, 0.12, 0.09, 0.48, 0.28
  )),
  sample_id = rep(c("SAMPLE_0008", "SAMPLE_0054"), each = 15),
  replicate = rep(c("A1", "A1"), each = 15)
)

# Create toy data for H3K27ac (enhancer/promoter marker)
h3k27ac_enrichment <- tibble(
  Chromatin_State = rep(chromatin_states, 2),
  sample_id_rep = rep(c("Proj1_A1_H3K27ac_1_SAMPLE_0008", 
                        "Proj1_A1_H3K27ac_1_SAMPLE_0054"), each = 15),
  mean_rpgc_val = c(
    # Sample 0008 - high at promoters and enhancers
    2.0, 1.5, 1.0, 0.9, 0.6, 0.5, 1.8, 1.6, 2.1, 1.9, 1.2, 0.2, 0.1, 0.4, 1.0,
    # Sample 0054 - similar pattern
    1.9, 1.4, 0.95, 0.85, 0.55, 0.45, 1.7, 1.5, 2.0, 1.8, 1.1, 0.18, 0.12, 0.38, 0.95
  ),
  mean_rpgc_text = sprintf("%.2f", c(
    2.0, 1.5, 1.0, 0.9, 0.6, 0.5, 1.8, 1.6, 2.1, 1.9, 1.2, 0.2, 0.1, 0.4, 1.0,
    1.9, 1.4, 0.95, 0.85, 0.55, 0.45, 1.7, 1.5, 2.0, 1.8, 1.1, 0.18, 0.12, 0.38, 0.95
  )),
  sample_id = rep(c("SAMPLE_0008", "SAMPLE_0054"), each = 15),
  replicate = rep(c("A1", "A1"), each = 15)
)

# Create toy data for 5mC (DNA methylation)
methylation_enrichment <- tibble(
  Chromatin_State = rep(chromatin_states, 2),
  sample_id_rep = rep(c("Proj1_A1_5mC_1_SAMPLE_0008", 
                        "Proj1_A1_5mC_1_SAMPLE_0054"), each = 15),
  mean_rpgc_val = c(
    # Sample 0008 - generally high, lower at promoters
    0.8, 0.5, 0.6, 0.7, 1.2, 1.1, 0.9, 0.9, 0.8, 0.8, 0.9, 1.0, 1.5, 0.4, 0.7,
    # Sample 0054 - similar pattern
    0.75, 0.48, 0.58, 0.68, 1.15, 1.08, 0.88, 0.88, 0.78, 0.78, 0.88, 0.98, 1.48, 0.38, 0.68
  ),
  mean_rpgc_text = sprintf("%.2f", c(
    0.8, 0.5, 0.6, 0.7, 1.2, 1.1, 0.9, 0.9, 0.8, 0.8, 0.9, 1.0, 1.5, 0.4, 0.7,
    0.75, 0.48, 0.58, 0.68, 1.15, 1.08, 0.88, 0.88, 0.78, 0.78, 0.88, 0.98, 1.48, 0.38, 0.68
  )),
  sample_id = rep(c("SAMPLE_0008", "SAMPLE_0054"), each = 15),
  replicate = rep(c("A1", "A1"), each = 15)
)

# Combine into named list
enrichment_results <- list(
  H3K4me3 = h3k4me3_enrichment,
  H3K27ac = h3k27ac_enrichment,
  `5mC` = methylation_enrichment
)

# ---- Create toy profile_results ----
# Structure: nested list [marker][feature_type] containing profile data

# Generate example profile data around TSS
# Index represents distance from feature center (e.g., -2000 to +2000 bp)
n_bins <- 100
index_vals <- seq(-2000, 2000, length.out = n_bins)

# Helper function to create profile data
create_profile_data <- function(marker, sample_ids, peak_height, peak_width) {
  n_samples <- length(sample_ids)
  
  # Create profile shape: Gaussian-like peak at center
  profile_shape <- peak_height * exp(-(index_vals^2) / (2 * peak_width^2))
  
  # Add some variation between samples and noise
  profiles <- lapply(seq_along(sample_ids), function(i) {
    variation <- rnorm(1, mean = 1, sd = 0.1)
    noise <- rnorm(n_bins, mean = 0, sd = peak_height * 0.05)
    profile_shape * variation + noise
  })
  
  # Flatten to long format
  tibble(
    mean_rpgc_val = unlist(profiles),
    sderror_rpgc_val = abs(rnorm(n_bins * n_samples, mean = peak_height * 0.1, sd = peak_height * 0.02)),
    median = unlist(lapply(profiles, function(p) p * 0.95)),
    index = rep(index_vals, n_samples),
    sample = rep(sample_ids, each = n_bins),
    min_error = NA_real_,
    max_error = NA_real_,
    marker = marker,
    sample_id = rep(gsub(".*_(SAMPLE_[0-9]+).*", "\\1", sample_ids), each = n_bins),
    sample_id_batch = rep(sample_ids, each = n_bins)
  )
}

# H3K4me3 - sharp peak at TSS (protein coding genes)
h3k4me3_protein_coding <- create_profile_data(
  marker = "H3K4me3",
  sample_ids = c("Proj1_A1_H3K4me3_1_SAMPLE_0008", 
                 "Proj1_A1_H3K4me3_1_SAMPLE_0054"),
  peak_height = 2.5,
  peak_width = 400
)

# H3K27ac - broader peak at TSS/enhancers
h3k27ac_protein_coding <- create_profile_data(
  marker = "H3K27ac",
  sample_ids = c("Proj1_A1_H3K27ac_1_SAMPLE_0008", 
                 "Proj1_A1_H3K27ac_1_SAMPLE_0054"),
  peak_height = 2.0,
  peak_width = 600
)

# 5mC - depletion at CpG islands
methylation_cpg <- create_profile_data(
  marker = "5mC",
  sample_ids = c("Proj1_A1_5mC_1_SAMPLE_0008", 
                 "Proj1_A1_5mC_1_SAMPLE_0054"),
  peak_height = -0.5,  # Negative = depletion
  peak_width = 800
) %>%
  mutate(
    # Shift baseline up and invert (high background, low at islands)
    mean_rpgc_val = 1.2 - mean_rpgc_val,
    median = 1.2 - median
  )

# Combine into nested list structure
profile_results <- list(
  H3K4me3 = list(
    protein_coding = h3k4me3_protein_coding
  ),
  H3K27ac = list(
    protein_coding = h3k27ac_protein_coding
  ),
  `5mC` = list(
    cpg_islands = methylation_cpg
  )
)

# ---- Save as R data objects ----
usethis::use_data(enrichment_results, overwrite = TRUE)
usethis::use_data(profile_results, overwrite = TRUE)

message("\n✓ Toy enrichment data objects created successfully!")
message("Objects saved to data/ directory:")
message("  - enrichment_results (", length(enrichment_results), " markers)")
message("  - profile_results (", length(profile_results), " markers)")
message("\nUsers can load with:")
message("  data(enrichment_results)")
message("  data(profile_results)")
