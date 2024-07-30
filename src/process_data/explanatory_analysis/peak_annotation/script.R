
library(IRanges)
library(GenomicRanges)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(stringr)
library(anndata)
library(tibble)


## VIASH START
par <- list(
  multiomics_atac = "resources/grn-benchmark/multiomics_atac.h5ad",
  annot_peak_database = "resources/grn-benchmark/supp/annot_peak_database.csv"
)
print(par)

adata_atac <- read_h5ad(par$multiomics_atac)

format_peak <- function(peaks) {
  formatted_peaks <- vector("character", length(peaks))
  for (i in seq_along(peaks)) {
    parts <- unlist(strsplit(peaks[i], "[:\\-_]"))
    chr <- parts[1]
    peak <- parts[2]
    formatted_peaks[i] <- paste0(chr, ":", peak)
  }
  return(formatted_peaks)
}

# Get the peak names and format them
atac_peaks <- format_peak(rownames(adata_atac$var))

# Create a data frame with chromosome and range information
peaks <- tibble(
  chr = sapply(atac_peaks, function(x) strsplit(x, ":")[[1]][1]),
  range = sapply(atac_peaks, function(x) paste(strsplit(x, ":")[[1]][2], collapse = ""))
)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks = GRanges(peaks$chr, IRanges(peaks$range))
peaks_annotated = suppressMessages(ChIPseeker::annotatePeak(
    peaks,
    tssRegion = c(-1000, 1000), # extended from -5kb to 5
    TxDb = txdb,
    level = "transcript", 
    assignGenomicAnnotation = TRUE,  # the default
    genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                "Downstream", "Intergenic"),  # the default
    annoDb = NULL,
    sameStrand = FALSE, # the default
    ignoreOverlap = FALSE, # the default
    ignoreUpstream = FALSE, # the default
    ignoreDownstream = FALSE, # the default
    overlap = "TSS", # the default
    verbose = TRUE # the default
))
peaks_annotated_df = as.data.frame(peaks_annotated)

# Define the mapping dictionary
map_ <- c('Intron' = 'Intron', 
          'Exon' = 'Exon', 
          'Promoter' = 'Promoter', 
          'Distal' = 'Distal Intergenic', 
          "3'" = "3' UTR", 
          'Downstream' = 'Downstream (<=300)', 
          "5'" = "5' UTR")

# Split the annotation column and map the values
peaks_annotated_df$ann <- sapply(str_split(peaks_annotated_df$annotation, ' ', simplify = TRUE)[, 1], function(x) map_[x])

# # Create the peaks string
peaks_annotated_df$peaks <- paste0(peaks_annotated_df$seqnames, ':', peaks_annotated_df$start, '-', peaks_annotated_df$end)

# # Create a new data frame with the required columns
df <- data.frame(annotation = peaks_annotated_df$ann, peak = peaks_annotated_df$peaks, stringsAsFactors = FALSE)

write.csv(df, file = par$annot_peak_database , row.names = TRUE)

