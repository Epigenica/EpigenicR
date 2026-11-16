#' Build non-overlapping genomic feature sets from a GTF
#'
#' Reads a gene annotation (GTF) and splits the genome into feature sets that
#' do not overlap each other: TSS, promoter, 5' UTR, 3' UTR, exon, intron and
#' intergenic. The result is a named list of \code{GRanges} that can be passed
#' as \code{annotations} to \code{\link{create_epk}} or
#' \code{\link{add_features_to_epk}}, which read the BigWig signal over each
#' region. The resulting experiments can then be plotted with
#' \code{\link{plot_features_stacked_bar}}.
#'
#' @param gtf_file Path to a GTF file (GENCODE or Ensembl style; may be
#'   gzipped), e.g. the file downloaded by \code{\link{ensure_gtf_and_beds}}.
#' @param chrom_sizes Chromosome lengths, needed to define intergenic regions.
#'   One of: a named numeric vector (\code{c(chr1 = 248956422, ...)}), the path
#'   to a BigWig file (lengths are read from its header), or the path to a
#'   two-column \code{chrom.sizes} file. If \code{NULL} (default), each
#'   chromosome is assumed to end at its last annotated gene, so the
#'   intergenic set misses the chromosome ends.
#' @param tss_flank Number of bases on each side of the transcription start
#'   site that form the \code{tss} set. Default \code{500}.
#' @param promoter_upstream Number of bases upstream of the TSS that form the
#'   \code{promoter} set (the part not already in \code{tss}). Default
#'   \code{2000}. Use \code{0} to leave out the promoter set.
#' @param protein_coding_only Logical; use only protein-coding genes. Regions
#'   of other genes then count as intergenic. Default \code{TRUE}.
#' @param chr_pattern Regular expression selecting the chromosomes to keep.
#'   Default keeps the autosomes, X and Y (with or without the \code{chr}
#'   prefix).
#'
#' @return A named list of \code{GRanges}, in priority order: \code{tss},
#'   \code{promoter}, \code{utr5}, \code{utr3}, \code{exon}, \code{intron},
#'   \code{intergenic}. Each has a \code{feature_id} metadata column. Empty
#'   sets are dropped.
#'
#' @details
#' A base that fits several definitions is assigned to the first matching set
#' in the priority order above, so every base of the kept chromosomes belongs
#' to exactly one set. Strand is ignored when resolving overlaps. For example,
#' \code{exon} holds only exonic bases that are not within a TSS window,
#' promoter or UTR, and \code{intron} holds gene-body bases not covered by any
#' exon of any transcript.
#'
#' GENCODE GTFs label both UTR types as \code{UTR}; these are split into 5'
#' and 3' by their position relative to the transcript's coding region.
#'
#' @examples
#' \dontrun{
#' sets <- build_genomic_feature_sets(
#'   gtf_file    = "data/gencode.v38.annotation.gtf",
#'   chrom_sizes = "minute_output/bigwig/one_of_the_files.bw"
#' )
#' epk <- add_features_to_epk(
#'   epk, annotations = sets,
#'   pipeline_output_path = "/path/to/project",
#'   bigwig_scale = "scaled", replicate_mode = "replicates"
#' )
#' plot_features_stacked_bar(epk, features = names(sets), metric = "signal")
#' }
#'
#' @export
build_genomic_feature_sets <- function(gtf_file,
                                       chrom_sizes = NULL,
                                       tss_flank = 500,
                                       promoter_upstream = 2000,
                                       protein_coding_only = TRUE,
                                       chr_pattern = "^(chr)?([0-9]+|X|Y)$") {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required. Install it with BiocManager::install('rtracklayer').")
  }
  if (!is.character(gtf_file) || length(gtf_file) != 1 || !file.exists(gtf_file)) {
    stop("'gtf_file' must be the path to an existing GTF file.")
  }
  if (tss_flank < 0 || promoter_upstream < 0) {
    stop("'tss_flank' and 'promoter_upstream' must be 0 or greater.")
  }

  gtf <- rtracklayer::import(gtf_file, format = "gtf")
  gtf <- gtf[grepl(chr_pattern, as.character(GenomicRanges::seqnames(gtf)))]
  if (length(gtf) == 0) {
    stop("No GTF records left after applying 'chr_pattern'.")
  }

  meta <- S4Vectors::mcols(gtf)
  if (protein_coding_only) {
    type_col <- intersect(c("gene_type", "gene_biotype"), names(meta))
    if (length(type_col) == 0) {
      stop(
        "The GTF has no 'gene_type' or 'gene_biotype' attribute. ",
        "Set protein_coding_only = FALSE."
      )
    }
    gtf <- gtf[!is.na(meta[[type_col[1]]]) & meta[[type_col[1]]] == "protein_coding"]
    meta <- S4Vectors::mcols(gtf)
    if (length(gtf) == 0) stop("No protein-coding records found in the GTF.")
  }

  rec_type <- as.character(meta$type)
  genes <- gtf[rec_type == "gene"]
  tx    <- gtf[rec_type == "transcript"]
  exons <- gtf[rec_type == "exon"]
  if (length(tx) == 0) tx <- genes
  if (length(genes) == 0 || length(exons) == 0) {
    stop("The GTF must contain 'gene' and 'exon' records.")
  }

  utrs <- .split_utrs(gtf)

  # Strand-aware windows around the transcription start site
  tss_gr <- GenomicRanges::promoters(tx, upstream = tss_flank, downstream = tss_flank + 1)
  prom_gr <- if (promoter_upstream > 0) {
    GenomicRanges::promoters(tx, upstream = promoter_upstream, downstream = 0)
  } else {
    tx[0]
  }

  # Chromosomes, in natural order (1..22, X, Y)
  chroms <- unique(as.character(GenomicRanges::seqnames(genes)))
  chr_key <- sub("^chr", "", chroms)
  chr_rank <- suppressWarnings(as.integer(chr_key))
  chr_rank[is.na(chr_rank)] <- 1000L + match(chr_key[is.na(chr_rank)], c("X", "Y"))
  chroms <- chroms[order(chr_rank)]

  genome_gr <- .genome_ranges(chrom_sizes, chroms, list(genes, tss_gr, prom_gr))

  candidates <- list(
    tss      = tss_gr,
    promoter = prom_gr,
    utr5     = utrs$utr5,
    utr3     = utrs$utr3,
    exon     = exons,
    intron   = genes,
    intergenic = genome_gr
  )

  # Assign each base to the first set that claims it
  claimed <- .plain_ranges(genome_gr[0], chroms)
  sets <- list()
  for (nm in names(candidates)) {
    piece <- .plain_ranges(candidates[[nm]], chroms)
    piece <- GenomicRanges::intersect(piece, genome_gr)
    piece <- GenomicRanges::setdiff(piece, claimed)
    claimed <- GenomicRanges::reduce(c(claimed, piece))
    if (length(piece) > 0) {
      S4Vectors::mcols(piece)$feature_id <- paste0(nm, "_", seq_along(piece))
      sets[[nm]] <- piece
    }
  }

  message(
    "Feature sets built (regions; % of genome): ",
    paste(
      vapply(names(sets), function(nm) {
        sprintf(
          "%s %d; %.1f%%", nm, length(sets[[nm]]),
          100 * sum(as.numeric(GenomicRanges::width(sets[[nm]]))) /
            sum(as.numeric(GenomicRanges::width(genome_gr)))
        )
      }, character(1)),
      collapse = " | "
    )
  )

  sets
}

#' Split UTR records into 5' and 3' UTRs
#' @keywords internal
.split_utrs <- function(gtf) {
  rec_type <- as.character(S4Vectors::mcols(gtf)$type)

  # Ensembl-style GTFs label the two UTR types directly
  if (any(rec_type %in% c("five_prime_utr", "three_prime_utr"))) {
    return(list(
      utr5 = gtf[rec_type == "five_prime_utr"],
      utr3 = gtf[rec_type == "three_prime_utr"]
    ))
  }

  utr <- gtf[rec_type == "UTR"]
  cds <- gtf[rec_type == "CDS"]
  if (length(utr) == 0 || length(cds) == 0) {
    return(list(utr5 = gtf[0], utr3 = gtf[0]))
  }

  # GENCODE-style: decide by position relative to the transcript's coding region
  cds_tx <- as.character(S4Vectors::mcols(cds)$transcript_id)
  cds_start <- tapply(GenomicRanges::start(cds), cds_tx, min)
  cds_end   <- tapply(GenomicRanges::end(cds), cds_tx, max)

  utr_tx <- as.character(S4Vectors::mcols(utr)$transcript_id)
  before_cds <- GenomicRanges::end(utr) < cds_start[utr_tx]
  after_cds  <- GenomicRanges::start(utr) > cds_end[utr_tx]
  minus <- as.character(GenomicRanges::strand(utr)) == "-"

  is_utr5 <- ifelse(minus, after_cds, before_cds)
  is_utr3 <- ifelse(minus, before_cds, after_cds)

  list(
    utr5 = utr[!is.na(is_utr5) & is_utr5],
    utr3 = utr[!is.na(is_utr3) & is_utr3]
  )
}

#' Strip strand and metadata, keep the chosen chromosomes, merge overlaps
#' @keywords internal
.plain_ranges <- function(gr, chroms) {
  chr <- as.character(GenomicRanges::seqnames(gr))
  keep <- chr %in% chroms & GenomicRanges::end(gr) >= 1
  gr <- gr[keep]
  out <- GenomicRanges::GRanges(
    seqnames = factor(chr[keep], levels = chroms),
    ranges = IRanges::IRanges(
      start = pmax(GenomicRanges::start(gr), 1L),
      end = GenomicRanges::end(gr)
    )
  )
  GenomicRanges::reduce(out)
}

#' Whole-chromosome ranges from chrom_sizes (or from the annotation extent)
#' @keywords internal
.genome_ranges <- function(chrom_sizes, chroms, annotation_list) {
  if (is.null(chrom_sizes)) {
    message(
      "'chrom_sizes' not given: each chromosome is taken to end at its last ",
      "annotated region, so the intergenic set misses the chromosome ends."
    )
    all_gr <- do.call(c, lapply(annotation_list, .plain_ranges, chroms = chroms))
    ends <- tapply(
      GenomicRanges::end(all_gr),
      as.character(GenomicRanges::seqnames(all_gr)),
      max
    )
    sizes <- stats::setNames(as.numeric(ends), names(ends))
  } else if (is.character(chrom_sizes) && length(chrom_sizes) == 1) {
    if (!file.exists(chrom_sizes)) {
      stop("'chrom_sizes' file not found: ", chrom_sizes)
    }
    if (grepl("\\.(bw|bigwig)$", chrom_sizes, ignore.case = TRUE)) {
      si <- rtracklayer::seqinfo(rtracklayer::BigWigFile(chrom_sizes))
      si_df <- as.data.frame(si)
      sizes <- stats::setNames(as.numeric(si_df$seqlengths), rownames(si_df))
    } else {
      tab <- utils::read.table(chrom_sizes, sep = "\t", stringsAsFactors = FALSE)
      sizes <- stats::setNames(as.numeric(tab[[2]]), tab[[1]])
    }
  } else if (is.numeric(chrom_sizes) && !is.null(names(chrom_sizes))) {
    sizes <- chrom_sizes
  } else {
    stop(
      "'chrom_sizes' must be NULL, a named numeric vector, a BigWig file path ",
      "or a chrom.sizes file path."
    )
  }

  missing_chr <- setdiff(chroms, names(sizes))
  if (length(missing_chr) > 0) {
    stop(
      "'chrom_sizes' has no length for chromosome(s): ",
      paste(missing_chr, collapse = ", ")
    )
  }

  GenomicRanges::GRanges(
    seqnames = factor(chroms, levels = chroms),
    ranges = IRanges::IRanges(start = 1, end = as.numeric(sizes[chroms]))
  )
}
