# Generate Example Plots for README
# This script creates example plot images from the toy dataset

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Source the plot function directly
source("R/functions.R")

# Load toy data from .rda files
load("data/toy_metadata.rda")
load("data/toy_stats_summary.rda")
load("data/toy_genes.rda")

# Assign to standard names
metadata <- toy_metadata
stats_summary <- toy_stats_summary

# Create output directory
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

# ---- 1. QC Plot - Multiple Statistics ----
cat("Generating QC plot with multiple statistics...\n")
qc_plot_multi <- plot_qc_stats(
  data = stats_summary,
  condition = "All",
  stats = c("final_mapped", "library_size", "percent_duplication"),
  marker_levels = unique(stats_summary$marker),
  engine = "ggplot",
  legend_position = "bottom",
  ncol = 3,
  sample_labeling = "sample_id"
)

ggsave(
  filename = "man/figures/qc_plot_example.png",
  plot = qc_plot_multi,
  width = 12,
  height = 4,
  dpi = 300,
  bg = "white"
)
cat("✓ Saved: man/figures/qc_plot_example.png\n")

# ---- 2. QC Plot - Single Statistic (Simpler) ----
cat("Generating simple QC plot...\n")
qc_plot_simple <- plot_qc_stats(
  data = stats_summary,
  condition = "All",
  stats = "final_mapped",
  marker_levels = unique(stats_summary$marker),
  engine = "ggplot",
  legend_position = "right",
  sample_labeling = "sample_id"
)

ggsave(
  filename = "man/figures/qc_plot_simple.png",
  plot = qc_plot_simple,
  width = 7,
  height = 4,
  dpi = 300,
  bg = "white"
)
cat("✓ Saved: man/figures/qc_plot_simple.png\n")

# ---- 3. Summary Statistics Table (as image) ----
cat("Generating summary statistics by marker...\n")

summary_stats <- stats_summary %>%
  group_by(marker) %>%
  summarise(
    n_samples = n(),
    mean_final_reads = round(mean(final_mapped) / 1e6, 2),
    mean_lib_size = round(mean(library_size) / 1e6, 2),
    mean_duplication = round(mean(percent_duplication), 1)
  ) %>%
  mutate(
    mean_final_reads = paste0(mean_final_reads, "M"),
    mean_lib_size = paste0(mean_lib_size, "M"),
    mean_duplication = paste0(mean_duplication, "%")
  )

# Print summary for reference
print(summary_stats)

cat("\n✓ All example plots generated successfully!\n")
cat("Images saved in: man/figures/\n")
