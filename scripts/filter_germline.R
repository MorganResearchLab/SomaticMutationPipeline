#! /usr/bin/env Rscript


suppressPackageStartupMessages({
library(dplyr)
library(VariantAnnotation)
library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(Homo.sapiens)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
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

# add gene symbols to this
txdb.hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

in.files <- list.files(opt$in_path, pattern=paste0(opt$sample, "_", opt$celltype), full.names=TRUE)
in.df.list <- list()

for(x in seq_along(in.files)){
  x.bed.df <- read.table(in.files[x], header=TRUE, comment.char="")
  x.bed.df <- x.bed.df[x.bed.df$Num_reads >= 2, , drop=FALSE]
  if(nrow(x.bed.df)){
    x.bed.df$ALT <- unlist(lapply(strsplit(x.bed.df$ALT_expected, split=",", fixed=TRUE), FUN=function(X) paste0(X[1])))
    x.file <- unlist(lapply(strsplit(in.files[x], split="/"), FUN=function(FX) paste0(FX[length(FX)])))
    x.samp <- gsub(x.file, pattern="(SRR[0-9]+)_(\\S*)(\\.single_cell_genotype\\.tsv)", replacement="\\1")
    x.ct <- gsub(x.file, pattern="(SRR[0-9]+)_(\\S*)_(chr[0-9]+)(\\.single_cell_genotype\\.tsv)", replacement="\\2")
    x.bed.df$CellType <- x.ct
    x.bed.df$Sample <- x.samp
    in.df.list[[in.files[x]]] <- x.bed.df
  }
}

bed.df <- do.call(rbind.data.frame, in.df.list)
rownames(bed.df) <- NULL
bed.df <- bed.df[, colnames(bed.df) %in% c("X.CHROM", "Start", "End", "REF", "ALT", "CellType", "Sample")]
bed.df <- bed.df %>% distinct(X.CHROM, Start, End, REF, ALT, .keep_all=TRUE)
mut.bed <- makeGRangesFromDataFrame(bed.df, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field="X.CHROM", start.field="Start",
                                    end.field="End")
mut.bed <- unique(mut.bed)


library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
# loop over somatic mutations and pull out which ones overlap with dbSNP entries.
seqlevelsStyle(mut.bed) <- "NCBI"
germline.snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, mut.bed)
germline.snps

not.germline <- mut.bed[!mut.bed %in% germline.snps]
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