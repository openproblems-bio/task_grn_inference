library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)


## VIASH START
par <- list(
  multiomics_rna = "resources_test/grn-benchmark/multiomics_r/rna.rds",
  multiomics_atac = "resources_test/grn-benchmark/multiomics_r/atac.rds",
  temp_dir =  "output/figr/",
  cell_topic = "resources_test/grn-benchmark/supp/cell_topic.csv",
  num_workers = 1,
  n_topics = 48,
  peak_gene = "output/figr/peak_gene.csv",
  prediction= "output/figr/prediction.csv"
)

print(par)
# meta <- list(
#   functionality_name = "my_method_r"
# )
## VIASH END
dir.create(par$temp_dir, recursive = TRUE, showWarnings = TRUE)

cellknn_func <- function(par) {
  ## load cell topic probabilities and create cell-cluster matrix
  cell_topic <- read.csv(paste0(par$cell_topic), row.names = 1)
  print(dim(cell_topic))
  # Derive cell kNN using this
  cellkNN <- get.knn(cell_topic, k=par$n_topics)$nn.index
  rownames(cellkNN) <- rownames(cell_topic)
  print(dim(cellkNN))
  print(paste0(par$temp_dir, "cellkNN.rds"))
  saveRDS(cellkNN, paste0(par$temp_dir, "cellkNN.rds"))
}
cellknn_func(par)

## Step1: Peak-gene association testing
peak_gene_func <- function(par){
  atac = readRDS(par$multiomics_atac)
  rna  = readRDS(par$multiomics_rna)
  cisCorr <- FigR::runGenePeakcorr(ATAC.se = atac,
                            RNAmat = rna,
                            genome = "hg38", # One of hg19, mm10 or hg38 
                            nCores = par$num_workers,
                            p.cut = NULL, # Set this to NULL and we can filter later
                            n_bg = 100)
  write.csv(cisCorr, paste0(par$temp_dir, "cisCorr.csv"), row.names = TRUE)
}

peak_gene_func(par)

## Step 2: create DORCs and smooth them 
dorc_genes_func <- function(par){
  cisCorr = read.csv(paste0(par$temp_dir, "cisCorr.csv"))
  cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

  atac = readRDS(par$multiomics_atac)
  rna  = readRDS(par$multiomics_rna)

  allGenes = unique(cisCorr.filt$Gene) 
  dorcMat <- getDORCScores(ATAC.se = atac, # Has to be same SE as used in previous step
                          dorcTab = cisCorr.filt,
                          geneList = allGenes,
                          nCores = par$num_workers)
  print(print(paste0(par$temp_dir, "cellkNN.rds")))
  cellkNN = readRDS(paste0(par$temp_dir, "cellkNN.rds"))
  # Smooth dorc scores using cell KNNs (k=n_topics)
  n_topics = par$n_topics
  common_cells <- intersect(rownames(cellkNN), colnames(rna))
  cellkNN = cellkNN[common_cells,]
  dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:n_topics], mat = dorcMat, nCores = 4)

  # Smooth RNA using cell KNNs
  # This takes longer since it's all genes
  RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:n_topics], mat = rna,nCores = 4)

  # get peak gene connection
  write.csv(cisCorr.filt, paste0(par$temp_dir, "cisCorr.filt.csv"))
  saveRDS(RNAmat.s, paste0(par$temp_dir, "RNAmat.s.RDS"))
  saveRDS(dorcMat.s, paste0(par$temp_dir, "dorcMat.s.RDS"))
}
dorc_genes_func(par)

## TF-gene associations
tf_gene_association_func <- function(par){
  cisCorr.filt = read.csv(paste0(par$temp_dir, "cisCorr.filt.csv"))
  RNAmat.s = readRDS(paste0(par$temp_dir, "RNAmat.s.RDS"))
  dorcMat.s = readRDS(paste0(par$temp_dir, "dorcMat.s.RDS"))

  atac = readRDS(par$multiomics_atac)
  figR.d <- runFigRGRN(ATAC.se = atac, # Must be the same input as used in runGenePeakcorr()
                      dorcTab = cisCorr.filt, # Filtered peak-gene associations
                      genome = "hg38",
                      dorcMat = dorcMat.s,
                      rnaMat = RNAmat.s, 
                      nCores = par$num_workers)

  write.csv(figR.d, paste0(par$temp_dir, "figR.d.csv"))
}
tf_gene_association_func(par)
extract_peak_gene_func <- function(par) {
  # Read the CSV file
  peak_gene_figr <- read.csv(file.path(par$temp_dir, "cisCorr.filt.csv"))
  
  # Calculate the number of peak ranges for each gene
  peak_gene_figr_n <- aggregate(PeakRanges ~ Gene, data = peak_gene_figr, length)
  
  # Calculate the max and median values
  max_peak_gene_figr_n <- max(peak_gene_figr_n$PeakRanges)
  median_peak_gene_figr_n <- median(peak_gene_figr_n$PeakRanges)
  
  # Print the results
  cat("In the peak-gene associations: number of CIS", length(unique(peak_gene_figr$PeakRanges)), 
      ", gene", length(unique(peak_gene_figr$Gene)), "\n")
  cat("Number of DORC genes", sum(peak_gene_figr_n$PeakRanges >= 10), "\n")
  
  # Select relevant columns and rename them
  peak_gene_figr <- peak_gene_figr[, c("PeakRanges", "Gene")]
  colnames(peak_gene_figr) <- c("peak", "target")
  
  # Write the result to a CSV file
  write.csv(peak_gene_figr, file = par$peak_gene, row.names = FALSE)
}
extract_peak_gene_func(par)

filter_figr_grn <- function(par) {
  # Read the CSV file
  figr_grn <- read.csv(file.path(par$temp_dir, "figR.d.csv"))
  
  # Filter based on enrichment
  figr_grn <- subset(figr_grn, Enrichment.P < 0.05)
  
  # Filter based on correlation
  figr_grn <- subset(figr_grn, Corr.P < 0.05)
  
  # Filter those that have a Score of 0
  figr_grn <- subset(figr_grn, Score != 0)
  
  # Subset columns
  figr_grn <- figr_grn[, c("Motif", "DORC", "Score")]
  
  # Reset row names (equivalent to resetting the index in Python)
  rownames(figr_grn) <- NULL
  
  # Rename columns
  colnames(figr_grn) <- c("source", "target", "weight")
  
  # Write the result to a CSV file
  write.csv(figr_grn, file = par$prediction, row.names = FALSE)
}


filter_figr_grn(par)