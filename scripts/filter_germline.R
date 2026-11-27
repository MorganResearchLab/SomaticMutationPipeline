#! /usr/bin/env Rscript


# Only the packages this script uses. The annotation packages it used to load
# (VariantAnnotation, GenomicFeatures, BSgenome.Hsapiens.UCSC.hg38, org.Hs.eg.db,
# Homo.sapiens, biomaRt, TxDb.Hsapiens.UCSC.hg38.knownGene, stringr) were never used
# and added a lot of memory per job. The dbSNP lookup below reads the SNPlocs data files
# directly, so SNPlocs/BSgenome are only loaded for the fallback.
suppressPackageStartupMessages({
library(dplyr)
library(GenomicRanges)
library(GenomeInfoDb)   # seqlevelsStyle()
library(Biostrings)
library(optparse)
})

parser <- OptionParser()
parser <- add_option(parser, c("--in_path"), type="character",
                     help="Directory containing genotype files")

parser <- add_option(parser, c("--sample"), type="character",
                     help="Sample ID to search on")

parser <- add_option(parser, c("--celltype"), type="character",
                     help="Cell type ID to search on")

parser <- add_option(parser, c("--output"), type="character",
                     help="Output file prefix of germline variants")

opt <- parse_args(parser)

# file latency problems on Maxwell so add a quick file inspection step.
system(paste0("ls -ls", opt$output))

# Match exactly <sample>_<celltype>_<chrom>.single_cell_genotype.tsv for this sample and
# cell type. An unanchored "<sample>_<celltype>" pattern also matched other cell types that
# share a prefix (e.g. NK -> NKT, NKhighresponse). Regex characters in the names are escaped,
# and the chromosome part allows any contig name (chr1-22, chrX, chrM, ...).
escape_regex <- function(s) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", s)
file.pattern <- paste0("^", escape_regex(opt$sample), "_", escape_regex(opt$celltype),
                       "_[^_]+\\.single_cell_genotype\\.tsv$")
in.files <- list.files(opt$in_path, pattern=file.pattern, full.names=TRUE)
if(length(in.files) == 0){
  stop(paste0("No genotype files matching ", file.pattern, " in ", opt$in_path))
}
in.df.list <- list()

for(x in seq_along(in.files)){
  x.bed.df <- read.table(in.files[x], header=TRUE, comment.char="")
  x.bed.df <- x.bed.df[x.bed.df$Num_reads >= 2, , drop=FALSE]
  if(nrow(x.bed.df)){
    x.bed.df$ALT <- unlist(lapply(strsplit(x.bed.df$ALT_expected, split=",", fixed=TRUE), FUN=function(X) paste0(X[1])))
    # Files are selected by exact sample and cell type above, so take both from the
    # arguments (the old filename regexes assumed SRR IDs and numeric chromosomes).
    x.bed.df$CellType <- opt$celltype
    x.bed.df$Sample <- opt$sample
    in.df.list[[in.files[x]]] <- x.bed.df
  }
}

# Write a header-only output and stop without error. Used when a cell type has no
# genotype rows with >= 2 reads, or when every candidate overlaps dbSNP (previously
# the script crashed in dplyr::distinct / downstream and failed the Snakemake job).
write_empty_and_quit <- function(reason) {
  message(reason, " - writing header-only output to ", opt$output)
  empty <- data.frame(seqnames=character(), start=integer(), end=integer(), width=integer(),
                      strand=character(), REF=character(), ALT=character(), CellType=character(),
                      Sample=character(), original_ref=character(), original_alt=character(),
                      was_flipped=logical(), Mutation=character())
  write.table(empty, file=opt$output, quote=FALSE, row.names=FALSE, sep="\t")
  quit(save="no", status=0)
}

if(length(in.df.list) == 0){
  write_empty_and_quit(paste0("No genotype rows with Num_reads >= 2 for ", opt$sample, " ", opt$celltype))
}

bed.df <- do.call(rbind.data.frame, in.df.list)
rownames(bed.df) <- NULL
bed.df <- bed.df[, colnames(bed.df) %in% c("X.CHROM", "Start", "End", "REF", "ALT", "CellType", "Sample")]
bed.df <- bed.df %>% distinct(X.CHROM, Start, End, REF, ALT, .keep_all=TRUE)
mut.bed <- makeGRangesFromDataFrame(bed.df, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field="X.CHROM", start.field="Start",
                                    end.field="End")
mut.bed <- unique(mut.bed)


# Flag candidate sites that are dbSNP155 positions.
#
# snpsByOverlaps() loads every dbSNP position on each chromosome it touches, even for one
# query site (~3.8 GB peak for chr22; chr1 is ~5x larger). The SNPlocs package stores
# dbSNP as an OnDiskLongTable: 4760 position-sorted batches of 200,000 SNPs, with
# spatial_index.rds giving each batch's genomic range and c001/bNNNNN.rds its positions.
# Reading only the batches that cover the candidate sites gives the same answer with a
# few MB of memory. If the package layout is not as expected, fall back to snpsByOverlaps().
#
# A site counts as germline when it is a single base whose position is in dbSNP, which is
# what the previous `mut.bed %in% snpsByOverlaps(...)` test did (exact range match against
# 1-bp SNP positions, any allele).
dbsnp_position_hits <- function(gr) {
  d <- system.file("extdata", package="SNPlocs.Hsapiens.dbSNP155.GRCh38")
  si.file <- file.path(d, "spatial_index.rds")
  if(!nzchar(d) || !file.exists(si.file)) return(NULL)
  si <- readRDS(si.file)
  batch.file <- function(b) file.path(d, "c001", sprintf("b%05d.rds", b))
  if(!file.exists(batch.file(length(si)))) return(NULL)
  is.hit <- logical(length(gr))
  ov <- findOverlaps(gr, si, ignore.strand=TRUE)
  for(b in unique(subjectHits(ov))){
    q <- unique(queryHits(ov)[subjectHits(ov) == b])
    pos <- readRDS(batch.file(b))
    is.hit[q] <- is.hit[q] | (start(gr)[q] %in% pos)
  }
  is.hit & width(gr) == 1L
}

seqlevelsStyle(mut.bed) <- "NCBI"
is.germline <- dbsnp_position_hits(mut.bed)
if(is.null(is.germline)){
  message("SNPlocs on-disk layout not recognised - falling back to snpsByOverlaps() (high memory)")
  suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP155.GRCh38))
  germline.snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, mut.bed)
  is.germline <- mut.bed %in% germline.snps
}
message(sum(is.germline), " of ", length(mut.bed), " candidate sites are dbSNP155 positions")

not.germline <- mut.bed[!is.germline]
if(length(not.germline) == 0){
  write_empty_and_quit(paste0("All ", length(mut.bed), " candidate sites overlap dbSNP for ", opt$sample, " ", opt$celltype))
}
seqlevelsStyle(not.germline) <- "UCSC"

check_strand_consistency <- function(gr) {
  # Get the reference and alternate alleles
  ref_alleles <- gr$REF
  alt_alleles <- gr$ALT
  
  # Convert to DNAString objects
  ref_dna <- DNAStringSet(ref_alleles)
  alt_dna <- DNAStringSet(alt_alleles)
  
  # Function to check if sequence pairs are reverse complements
  is_reverse_complement <- function(seq1, seq2) {
    rc_seq1 <- reverseComplement(seq1)
    return(as.character(rc_seq1) == as.character(seq2))
  }
  
  # Check each REF/ALT pair
  inconsistent_pairs <- mapply(function(ref, alt) {
    # Check if either the direct sequence or its reverse complement matches
    direct_match <- as.character(ref) == as.character(alt)
    rc_match <- is_reverse_complement(ref, alt)
    
    # If neither matches, it's consistent (different alleles on same strand)
    # If either matches, it might indicate a strand issue
    return(direct_match || rc_match)
  }, ref_dna, alt_dna)
  
  # Add results to the GRanges object
  gr$possible_strand_issue <- inconsistent_pairs
  
  # Create a summary
  summary <- list(
    total_variants = length(gr),
    possible_strand_issues = sum(inconsistent_pairs),
    consistent_variants = sum(!inconsistent_pairs)
  )
  
  return(list(
    granges = gr,
    summary = summary
  ))
}

# check the strand of all somatic mutations
result <- check_strand_consistency(not.germline)


standardize_variant_strands <- function(gr, target_strand = '+') {
  # Validate input
  if (!all(c('REF', 'ALT') %in% names(mcols(gr)))) {
    stop("GRanges object must contain REF and ALT columns")
  }
  if (!(target_strand %in% c('+', '-'))) {
    stop("target_strand must be either '+' or '-'")
  }
  
  # Create DNAStringSet objects for manipulation
  ref_dna <- DNAStringSet(gr$REF)
  alt_dna <- DNAStringSet(gr$ALT)
  
  # Function to check if sequence needs to be reverse complemented
  needs_rc <- function(ref, alt) {
    rc_ref <- reverseComplement(ref)
    rc_alt <- reverseComplement(alt)
    
    # Check if reverse complement version makes more sense
    # (i.e., if REF/ALT are reverse complements of each other)
    direct_pair <- c(as.character(ref), as.character(alt))
    rc_pair <- c(as.character(rc_ref), as.character(rc_alt))
    
    # Return TRUE if the reverse complement version should be used
    return(paste(sort(rc_pair), collapse="") < paste(sort(direct_pair), collapse=""))
  }
  
  # Identify variants that need to be reverse complemented
  to_flip <- mapply(needs_rc, ref_dna, alt_dna)
  
  # Create new GRanges object
  gr_new <- gr
  
  # Update sequences that need to be flipped
  if (any(to_flip)) {
    gr_new$REF[to_flip] <- as.character(reverseComplement(ref_dna[to_flip]))
    gr_new$ALT[to_flip] <- as.character(reverseComplement(alt_dna[to_flip]))
  }
  
  # Set all strands to target strand
  strand(gr_new) <- target_strand
  
  # Add metadata about changes
  gr_new$original_ref <- gr$REF
  gr_new$original_alt <- gr$ALT
  gr_new$was_flipped <- to_flip
  
  # Create summary
  summary <- list(
    total_variants = length(gr),
    variants_flipped = sum(to_flip),
    variants_unchanged = sum(!to_flip)
  )
  
  return(list(
    granges = gr_new,
    summary = summary
  ))
}

# Function to validate the standardization
validate_standardization <- function(gr_original, gr_standardized) {
  # Check that all variants are now on the same strand
  all_same_strand <- all(strand(gr_standardized) == strand(gr_standardized)[1])
  
  # Check that REF/ALT pairs are consistent
  ref_dna_orig <- DNAStringSet(gr_original$REF)
  alt_dna_orig <- DNAStringSet(gr_original$ALT)
  ref_dna_std <- DNAStringSet(gr_standardized$REF)
  alt_dna_std <- DNAStringSet(gr_standardized$ALT)
  
  # For each variant, check that the standardized version represents 
  # the same change as the original
  consistent_changes <- mapply(function(ref_o, alt_o, ref_s, alt_s) {
    orig_pair <- sort(c(as.character(ref_o), as.character(alt_o)))
    std_pair <- sort(c(as.character(ref_s), as.character(alt_s)))
    rc_std_pair <- sort(c(
      as.character(reverseComplement(DNAString(ref_s))),
      as.character(reverseComplement(DNAString(alt_s)))
    ))
    
    return(
      identical(orig_pair, std_pair) ||
      identical(orig_pair, rc_std_pair)
    )
  }, ref_dna_orig, alt_dna_orig, ref_dna_std, alt_dna_std)
  
  return(list(
    all_same_strand = all_same_strand,
    all_changes_consistent = all(consistent_changes),
    inconsistent_variants = which(!consistent_changes)
  ))
}


# Standardize all variants to the positive strand
pos.strand.grt <- standardize_variant_strands(not.germline, target_strand = '+')

# Print summary of changes
print(pos.strand.grt$summary)

# Validate the standardization
pos.strand.validation <- validate_standardization(not.germline, pos.strand.grt$granges)
print(pos.strand.validation)

# Access the standardized GRanges object
standardized.muts <- pos.strand.grt$granges

# Check which variants were flipped
flipped_variants <- standardized.muts[standardized.muts$was_flipped]

mut.to.save <- as.data.frame(standardized.muts)
mut.to.save$Mutation <- paste(paste0(mut.to.save$seqnames, ":", paste(mut.to.save$start, mut.to.save$end, sep="-")),
                              mut.to.save$REF, mut.to.save$ALT, sep="_")
write.table(mut.to.save, file=opt$output,
            quote=FALSE, row.names=FALSE, sep="\t")
