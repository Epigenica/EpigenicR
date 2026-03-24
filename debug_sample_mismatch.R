# Debug: Check which samples each marker has

source("create_sample_metadata.R")

# Group by marker and show unique samples
cat("=== Samples per marker ===\n\n")
for (m in c("INPUT", "M2", "M3", "M4")) {
  marker_data <- sample_metadata[sample_metadata$marker == m, ]
  cat(sprintf("\n%s (%d files):\n", m, nrow(marker_data)))
  print(marker_data)
  cat("\nUnique sample_id values:\n")
  print(sort(unique(marker_data$sample_id)))
}

# Check which samples are present in all markers (excluding INPUT)
cat("\n\n=== Sample coverage across markers ===\n")
m2_samples <- unique(sample_metadata$sample_id[sample_metadata$marker == "M2"])
m3_samples <- unique(sample_metadata$sample_id[sample_metadata$marker == "M3"])
m4_samples <- unique(sample_metadata$sample_id[sample_metadata$marker == "M4"])

cat("\nM2 samples:", length(m2_samples), "\n")
print(sort(m2_samples))

cat("\nM3 samples:", length(m3_samples), "\n")
print(sort(m3_samples))

cat("\nM4 samples:", length(m4_samples), "\n")
print(sort(m4_samples))

# Find common samples across all markers
common_samples <- Reduce(intersect, list(m2_samples, m3_samples, m4_samples))
cat("\n\nCommon samples across M2, M3, M4:", length(common_samples), "\n")
print(sort(common_samples))

# Find missing samples
cat("\n\n=== Missing samples ===\n")
all_samples <- unique(c(m2_samples, m3_samples, m4_samples))
for (sample in sort(all_samples)) {
  present <- c(
    if (sample %in% m2_samples) "M2" else "-",
    if (sample %in% m3_samples) "M3" else "-",
    if (sample %in% m4_samples) "M4" else "-"
  )
  cat(sprintf("%s: %s\n", sample, paste(present, collapse = " ")))
}
