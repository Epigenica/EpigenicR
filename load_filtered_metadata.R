# Complete working script: creates bw_files_filtered and sample_metadata_filtered
# for use with create_epk()

# Load original data
source("create_sample_metadata.R")

# Identify common samples across M2, M3, M4
m2_samples <- unique(sample_metadata$sample_id[sample_metadata$marker == "M2"])
m3_samples <- unique(sample_metadata$sample_id[sample_metadata$marker == "M3"])
m4_samples <- unique(sample_metadata$sample_id[sample_metadata$marker == "M4"])

common_samples <- Reduce(intersect, list(m2_samples, m3_samples, m4_samples))

# Filter to only:
# 1. Common samples across all markers
# 2. Exclude INPUT (will be excluded by markers_to_exclude anyway)
# 3. Keep only unscaled files (drop scaled duplicates)
keep_idx <- (sample_metadata$sample_id %in% common_samples) & 
            (sample_metadata$marker != "INPUT") &
            (!grepl("\\.scaled\\.bw$", bw_files))

bw_files_filtered <- bw_files[keep_idx]
sample_metadata_filtered <- sample_metadata[keep_idx, ]

cat("✓ Created bw_files_filtered:", length(bw_files_filtered), "files\n")
cat("✓ Created sample_metadata_filtered:", nrow(sample_metadata_filtered), "rows\n")
cat("✓ Markers:", paste(sort(unique(sample_metadata_filtered$marker)), collapse=", "), "\n")
cat("✓ Samples:", paste(sort(unique(sample_metadata_filtered$sample_id)), collapse=", "), "\n\n")

cat("Ready to use:\n")
cat("  bw_files_filtered\n")
cat("  sample_metadata_filtered\n")
