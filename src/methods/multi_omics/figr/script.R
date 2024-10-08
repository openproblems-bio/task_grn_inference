library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)


## VIASH START
par <- list(
  multiomics_rna_r = "resources/grn-benchmark/multiomics_rna.rds",
  multiomics_atac_r = "resources/grn-benchmark/multiomics_atac.rds",
  temp_dir =  "output/figr/",
  cell_topic = "resources/prior/cell_topic.csv",
  num_workers = 20,
  n_topics = 48,
  peak_gene = "output/figr/peak_gene.csv",
  prediction= "resources/grn_models/figr.csv"
)
print(par)
# meta <- list(
#   functionality_name = "my_method_r"
# )
## VIASH END
dir.create(par$temp_dir, recursive = TRUE, showWarnings = TRUE)

atac = readRDS(par$multiomics_atac_r)
colnames(atac) <- gsub("-", "", colnames(atac))

rna  = readRDS(par$multiomics_rna_r)
colnames(rna) <- gsub("-", "", colnames(rna))


cellknn_func <- function(par) {
  ## load cell topic probabilities and create cell-cluster matrix
  cell_topic <- read.csv(paste0(par$cell_topic), row.names = 1)
  print(dim(cell_topic))
  # Derive cell kNN using this
  cellkNN <- get.knn(cell_topic, k=par$n_topics)$nn.index
  rownames(cellkNN) <- rownames(cell_topic)
  saveRDS(cellkNN, paste0(par$temp_dir, "cellkNN.rds"))
}

## Step1: Peak-gene association testing
peak_gene_func <- function(par){
  
  common_cells <- intersect(colnames(atac), colnames(rna))
  rna = rna[,common_cells]
  atac = atac[,common_cells]

  cisCorr <- FigR::runGenePeakcorr(ATAC.se = atac,
                            RNAmat = rna,
                            genome = "hg38", # One of hg19, mm10 or hg38 
                            nCores = par$num_workers,
                            p.cut = NULL, # Set this to NULL and we can filter later
                            n_bg = 100)
  write.csv(cisCorr, paste0(par$temp_dir, "cisCorr.csv"), row.names = TRUE)
}

## Step 2: create DORCs and smooth them 
dorc_genes_func <- function(par){
  cisCorr = read.csv(paste0(par$temp_dir, "cisCorr.csv"))
  cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

  allGenes = unique(cisCorr.filt$Gene) 
  dorcMat <- getDORCScores(ATAC.se = atac, # Has to be same SE as used in previous step
                          dorcTab = cisCorr.filt,
                          geneList = allGenes,
                          nCores = par$num_workers)
  cellkNN = readRDS(paste0(par$temp_dir, "cellkNN.rds"))

  # Smooth dorc scores using cell KNNs (k=n_topics)
  n_topics = par$n_topics
  common_cells <- intersect(rownames(cellkNN), colnames(rna))

  cellkNN = cellkNN[common_cells,]
  dorcMat = dorcMat[,common_cells]
  cat('cellKNN dim:', dim(cellkNN), '\n')
  cat('dorcMat dim:', dim(dorcMat), '\n')
  cat('rna dim:', dim(rna), '\n')
  dorcMat.s <- smoothScoresNN(NNmat = cellkNN, mat = dorcMat, nCores = par$num_workers) 
  cat('dorcMat.s completed')
  # Smooth RNA using cell KNNs
  # This takes longer since it's all genes
  matching_indices <- match(colnames(rna), rownames(cellkNN))
  cellkNN_ordered <- cellkNN[matching_indices, ]


  RNAmat.s <- smoothScoresNN(NNmat = cellkNN_ordered, mat = rna, nCores = par$num_workers)
  # RNAmat.s <- rna
  cat('RNAmat.s completed')
  # get peak gene connection
  write.csv(cisCorr.filt, paste0(par$temp_dir, "cisCorr.filt.csv"))
  saveRDS(RNAmat.s, paste0(par$temp_dir, "RNAmat.s.RDS"))
  saveRDS(dorcMat.s, paste0(par$temp_dir, "dorcMat.s.RDS"))
}

## TF-gene associations
tf_gene_association_func <- function(par){
  cisCorr.filt = read.csv(paste0(par$temp_dir, "cisCorr.filt.csv"))
  RNAmat.s = readRDS(paste0(par$temp_dir, "RNAmat.s.RDS"))
  dorcMat.s = readRDS(paste0(par$temp_dir, "dorcMat.s.RDS"))

  figR.d <- runFigRGRN(ATAC.se = atac, # Must be the same input as used in runGenePeakcorr()
                      dorcTab = cisCorr.filt, # Filtered peak-gene associations
                      genome = "hg38",
                      dorcMat = dorcMat.s,
                      rnaMat = RNAmat.s, 
                      nCores = par$num_workers)

  write.csv(figR.d, paste0(par$temp_dir, "figR.d.csv"))
}

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

filter_figr_grn <- function(par) {
  # Read the CSV file
  figr_grn <- read.csv(file.path(par$temp_dir, "figR.d.csv"))

  # Filter those that have a Score of 0
  figr_grn <- subset(figr_grn, Score != 0)
  
  # Filter based on enrichment
  figr_grn <- subset(figr_grn, Enrichment.P < 0.05)
  
  # Filter based on correlation
  # figr_grn <- subset(figr_grn, Corr.P < 0.05)
  
  
  # Subset columns
  figr_grn <- figr_grn[, c("Motif", "DORC", "Score")]
  
  # Reset row names (equivalent to resetting the index in Python)
  rownames(figr_grn) <- NULL
  
  # Rename columns
  colnames(figr_grn) <- c("source", "target", "weight")
  
  # Write the result to a CSV file
  write.csv(figr_grn, file = par$prediction, row.names = FALSE)
}




cellknn_func(par)
print('1: cellknn_func finished')
peak_gene_func(par)
print('2: peak_gene_func finished')
dorc_genes_func(par)
print('3: dorc_genes_func finished')
tf_gene_association_func(par)
print('3: tf_gene_association_func finished')
extract_peak_gene_func(par)
print('4: extract_peak_gene_func finished')
filter_figr_grn(par)