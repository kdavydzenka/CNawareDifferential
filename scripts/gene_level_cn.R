setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/CRC/")

pkgs <- c("dplyr", "tidyr", "purr", "GenomicRanges", "IRanges", 
          "tibble", "biomaRt", "Homo.sapiens", "org.Hs.eg.db", "AnnotationHub")
sapply(pkgs, require, character.only = TRUE)

input_dir <- "seg_files"

seg_files <- list.files(input_dir, pattern = "\\.seg$", full.names = TRUE)

# Load GENCODE annotation 

ah <- AnnotationHub()
gc <- ah[["AH49556"]]  # Gencode v29 GFF3 for hg38
genes_gr <- keepStandardChromosomes(gc, pruning.mode = "coarse")

# Prepare mapping: gene coordinates and gene names
gene_ranges <- GenomicRanges::GRanges(seqnames = seqnames(genes_gr),
                                      ranges   = ranges(genes_gr),
                                      gene     = genes_gr$gene_name)

# Function to process a single SEG file -> gene-level CN

process_seg <- function(file, gene_ranges) {
  sample_id <- gsub(".*/|\\.seg$", "", file)
  
  seg <- read.delim(file, stringsAsFactors = FALSE)
  seg_gr <- GenomicRanges::GRanges(seqnames = seg$Chromosome,
                    ranges   = IRanges(seg$Start, seg$End),
                    score    = seg$Segment_Mean)
  
  # Overlap segments with genes
  hits <- IRanges::findOverlaps(seg_gr, gene_ranges)
  
  df <- tibble(
    GeneID = mcols(gene_ranges)$gene[subjectHits(hits)],
    Segment_Mean = mcols(seg_gr)$score[queryHits(hits)]
  ) %>%
    dplyr::group_by(GeneID) %>%
    dplyr::summarise("{sample_id}" := mean(Segment_Mean, na.rm = TRUE), .groups = "drop")
  
  return(df)
}

# Process all samples
gene_data_list <- lapply(seg_files, process_seg, gene_ranges = gene_ranges)

# Merge all into one gene-level matrix
merged_gene_data <- purr::reduce(gene_data_list, full_join, by = "GeneID")

saveRDS(merged_gene_data, file = "test_data/cn_gene_level.rds")