# create_metadata_df(bw_files = file.path(bigwig_dir, bw_df$bw_file))
smpl_mixed_name <- paste0(bw_df$sample_id, '_', bw_df$replicate)
bw_files <- file.path(bigwig_dir, bw_df$bw_file)
# length(bw_files)
# bw_df %>% filter(replicate == "pooled") %>% select(marker) %>% table()
# reading files per marker based the presence of the marker names within bigwig files and creating GRanges object to save ina list object and converting it to sparse matrix. 
assay_list_protein_coding <- list()
assay_list_cpg_islands <- list()
assay_list_tss_2k <- list()
assay_list_vista_enhancer <- list()
assay_list_E107_chromatin_state <- list()

excluded_markers <- c("INPUT")
markers_to_run <- setdiff(markers, excluded_markers)
markers_to_run[4] <- '5mC'
for(mk in markers_to_run){
  message("Processing marker: ", mk)
  idx <- grep(pattern = mk, x = bw_files, perl = T)
  bw_files_tmp <- bw_files[idx]
  idx <- grep(x = bw_files_tmp, pattern = '_pooled', perl = T, invert = T) # exclude pooled replicates
  bw_files_marker <- bw_files_tmp[idx]
  smpl_mixed_name_marker <- smpl_mixed_name[idx]
  # Protein_coding genes
  message("Processing protein coding genes for marker: ", mk)
  bw_gr_marker <- bw_loci(
    bwfiles = bw_files_marker,
    loci = genes_coord_protein_coding,
    labels=smpl_mixed_name_marker
  )
  colnames(mcols(bw_gr_marker)) <- gsub(x = smpl_mixed_name_marker, pattern = '\\.', replacement = '-')
  m <- as.matrix(S4Vectors::mcols(bw_gr_marker)[, smpl_mixed_name_marker, drop = FALSE])
  storage.mode(m) <- 'numeric'
  assay_list_protein_coding[[mk]] <- m

  # cpg islands
  message("Processing cpg islands for marker: ", mk)
  bw_gr_marker <- bw_loci(
    bwfiles = bw_files_marker,
    loci = cpg_islands,
    labels=smpl_mixed_name_marker
  )
  colnames(mcols(bw_gr_marker)) <- gsub(x = smpl_mixed_name_marker, pattern = '\\.', replacement = '-')
  m <- as.matrix(S4Vectors::mcols(bw_gr_marker)[, smpl_mixed_name_marker, drop = FALSE])
  storage.mode(m) <- 'numeric'
  assay_list_cpg_islands[[mk]] <- m

  # tss 2k
  message("Processing tss ±2kb for marker: ", mk)
  bw_gr_marker <- bw_loci(
    bwfiles = bw_files_marker,
    loci = tss_2kb_protein_coding,
    labels=smpl_mixed_name_marker
  )
  colnames(mcols(bw_gr_marker)) <- gsub(x = smpl_mixed_name_marker, pattern = '\\.', replacement = '-')
  m <- as.matrix(S4Vectors::mcols(bw_gr_marker)[, smpl_mixed_name_marker, drop = FALSE])
  storage.mode(m) <- 'numeric'
  assay_list_tss_2k[[mk]] <- m

  # vista enhancer
  bw_gr_marker <- bw_loci(
    bwfiles = bw_files_marker,
    loci = vista_enhancer,
    labels=smpl_mixed_name_marker
  )
  colnames(mcols(bw_gr_marker)) <- gsub(x = smpl_mixed_name_marker, pattern = '\\.', replacement = '-')
  m <- as.matrix(S4Vectors::mcols(bw_gr_marker)[, smpl_mixed_name_marker, drop = FALSE])
  storage.mode(m) <- 'numeric'
  assay_list_vista_enhancer[[mk]] <- m
  
  # E107 chromatin state
  message("Processing E107 chromatin state for marker: ", mk)
  bw_gr_marker <- bw_loci(
    bwfiles = bw_files_marker,
    loci = E107_chromatin_state,
    labels=smpl_mixed_name_marker
  )
  colnames(mcols(bw_gr_marker)) <- gsub(x = smpl_mixed_name_marker, pattern = '\\.', replacement = '-')
  m <- as.matrix(S4Vectors::mcols(bw_gr_marker)[, smpl_mixed_name_marker, drop = FALSE])
  storage.mode(m) <- 'numeric'
  assay_list_E107_chromatin_state[[mk]] <- m
}

# Creating SummarizedExperiment object
ref_cols <- smpl_mixed_name_marker
coldata <- S4Vectors::DataFrame(sample_id = ref_cols, row.names = ref_cols)

# Protein coding
se_protein_coding <- SummarizedExperiment::SummarizedExperiment(
  assays = assay_list_protein_coding,
  rowRanges = genes_coord_protein_coding,
  colData = coldata
)
rownames(se_protein_coding) <- genes_coord_protein_coding$gene_name

# cpg islands
se_cpg_islands<- SummarizedExperiment::SummarizedExperiment(
  assays = assay_list_cpg_islands,
  rowRanges = cpg_islands,
  colData = coldata
)
rownames(se_cpg_islands) <- cpg_islands$cpg_id

# tss 2kb
se_tss_2k<- SummarizedExperiment::SummarizedExperiment(
  assays = assay_list_tss_2k,
  rowRanges = tss_2kb_protein_coding,
  colData = coldata
)
rownames(se_tss_2k) <- tss_2kb_protein_coding$gene_name

# vista enhancer
se_vista_enhancer<- SummarizedExperiment::SummarizedExperiment(
  assays = assay_list_vista_enhancer,
  rowRanges = vista_enhancer,
  colData = coldata
)
rownames(se_vista_enhancer) <- vista_enhancer$enhancer

# E107 chromatin state
se_E107_chromatin_state<- SummarizedExperiment::SummarizedExperiment(
  assays = assay_list_E107_chromatin_state,
  rowRanges = E107_chromatin_state,
  colData = coldata
)
rownames(se_E107_chromatin_state) <- E107_chromatin_state$E107_chromatin_state

objlist <- list("protein_coding" = se_protein_coding, 
                "cpg" = se_cpg_islands, 
                "tss_2k" = se_tss_2k,
                "vista_enhancer" = se_vista_enhancer,
                "E107_chromatin_state" = se_E107_chromatin_state)

mse <- MultiAssayExperiment::MultiAssayExperiment(experiments = objlist, colData = coldata)

# Loading enrichment results paths
files.list_protein_coding <- list.files(file.path(proj_dir, "results/Enrichment", "protein_coding"), pattern = ".*_chromatin_state_dist\\.csv$", full.names = T,recursive = T)
files.list_cpg_islands <- file.path(proj_dir, "results/Enrichment", "CpG_islands", "5mC_chromatin_state_dist.csv")
files.list <- c(files.list_protein_coding, files.list_cpg_islands)
names(files.list) <- c(basename(dirname(files.list_protein_coding)), basename(dirname(files.list_cpg_islands)))
enrichment_results <- list(
  protein_coding = lapply(files.list_protein_coding, read.csv, sep = ",", header = T, stringsAsFactors = F),
  cpg_islands = lapply(files.list_cpg_islands, read.csv, sep = ",", header = T, stringsAsFactors = F)
)

files.list_protein_coding <- list.files(file.path(proj_dir, "results/Enrichment", "protein_coding"), pattern = ".*_profile_start_data\\.csv$", full.names = T,recursive = T)
files.list_cpg_islands <- file.path(proj_dir, "results/Enrichment", "CpG_islands", "methylation_profile_start_data.csv")
profile_results <- list(
  protein_coding = lapply(files.list_protein_coding, read.csv, sep = ",", header = T, stringsAsFactors = F),
  cpg_islands = lapply(files.list_cpg_islands, read.csv, sep = ",", header = T, stringsAsFactors = F)
)

epk <- structure(
  list(
    mse = mse,
    tables = list(
      stats_summary = stats_summary
    ),
    enrichment_results = list(
      chromatin_states = enrichment_results,
      enrichment_profile = profile_results
      ),
    
    provenance = list(created = Sys.time(),
                      session = sessionInfo())
  ),
  class = "EPK"
)
saveRDS( epk, file.path(results_qc, "project.epk.rds"))