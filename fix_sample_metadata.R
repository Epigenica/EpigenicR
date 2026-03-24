# Fix: Filter bw_files and sample_metadata to only include common samples across all markers

source("create_sample_metadata.R")

# Identify common samples across M2, M3, M4
m2_samples <- unique(sample_metadata$sample_id[sample_metadata$marker == "M2"])
m3_samples <- unique(sample_metadata$sample_id[sample_metadata$marker == "M3"])
m4_samples <- unique(sample_metadata$sample_id[sample_metadata$marker == "M4"])

common_samples <- Reduce(intersect, list(m2_samples, m3_samples, m4_samples))

cat("Common samples across all markers (M2, M3, M4):\n")
print(sort(common_samples))
cat("\n")

# Filter to only rows with common samples
# Also exclude INPUT markers and keep only unscaled files (not scaled)
keep_idx <- (sample_metadata$sample_id %in% common_samples) & 
            (sample_metadata$marker != "INPUT") &
            (!grepl("\\.scaled\\.bw$", bw_files))

# Create filtered versions
bw_files_filtered <- bw_files[keep_idx]
sample_metadata_filtered <- sample_metadata[keep_idx, ]

cat("Original: ", length(bw_files), "files\n")
cat("Filtered: ", length(bw_files_filtered), "files\n\n")

cat("Filtered sample counts by marker:\n")
print(table(sample_metadata_filtered$marker, sample_metadata_filtered$sample_id))

cat("\n\nFiltered metadata:\n")
print(sample_metadata_filtered)

# Export for use in create_epk
cat("\n\n=== Use these in your create_epk call ===\n")
cat("bw_files_filtered has", length(bw_files_filtered), "files\n")
cat("sample_metadata_filtered has", nrow(sample_metadata_filtered), "rows\n")
