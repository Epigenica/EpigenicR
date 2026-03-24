# Create sample_metadata data frame for the selected 29 BigWig files

bw_files <- c(
	"../../results/output_files_from_pipeline//bigwig/INP_Donor1_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/INP_Donor1_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/INP_Donor2_S1_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/INP_Donor2_S2_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/INP_Donor3_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/INP_Donor3_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/INP_Donor4_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M2_S2_Donor1_S1_rep2.hg19.scaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M2_S2_Donor1_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M2_S2_Donor1_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M2_S2_Donor2_S1_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M2_S2_Donor2_S2_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M2_S2_Donor3_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M2_S2_Donor3_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M2_S2_Donor4_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M2_S2_Donor4_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M3_S3_Donor1_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M3_S3_Donor2_S1_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M3_S3_Donor3_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M3_S3_Donor3_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M3_S3_Donor4_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M3_S3_Donor4_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M4_S4_Donor1_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M4_S4_Donor1_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M4_S4_Donor2_S1_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M4_S4_Donor2_S2_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M4_S4_Donor3_S1_rep2.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M4_S4_Donor3_S2_rep1.hg19.unscaled.bw",
	"../../results/output_files_from_pipeline//bigwig/M4_S4_Donor4_S2_rep1.hg19.unscaled.bw"
)

bw_base <- basename(bw_files)

sample_metadata <- data.frame(
	marker = ifelse(grepl("^INP_", bw_base), "INPUT", sub("^([^_]+)_.*$", "\\1", bw_base)),
	sample_id = sub(".*_(Donor[0-9]+_S[0-9]+)_(pooled|rep[0-9]+)\\.[^.]+\\.(scaled|unscaled)\\.bw$", "\\1", bw_base),
	replicate = sub(".*_(pooled|rep[0-9]+)\\.[^.]+\\.(scaled|unscaled)\\.bw$", "\\1", bw_base),
	stringsAsFactors = FALSE
)

print(sample_metadata)
cat("\nDimensions:", nrow(sample_metadata), "rows x", ncol(sample_metadata), "columns\n")
cat("\nCounts by marker:\n")
print(table(sample_metadata$marker))
