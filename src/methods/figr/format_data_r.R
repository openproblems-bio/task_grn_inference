library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(FigR)
library(optparse)


# Define the command line arguments
option_list <- list(
  make_option(c("--temp_dir"), type = "character", help = "Temporary directory to write the file"),
  make_option(c("--cell_topic"), type = "character", help = "Cell topic value")
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Extract all arguments
args <- commandArgs(trailingOnly = TRUE)
args_list <- list()
for (i in seq(1, length(args), by = 2)) {
  args_list[[gsub("--", "", args[i])]] <- args[i + 1]
}


## Load atac-seq and create summarizedexperiment
X <- readMM(paste0(par$temp_dir, "/scATAC/X_matrix.mtx"))
X <- t(X)
annotation_peak <- read.csv(paste0(par$temp_dir, "/scATAC/annotation_peak.csv"), row.names = 1)
annotation_cells <- read.csv(paste0(par$temp_dir, "/scATAC/annotation_cells.csv"), row.names = 1)

# Filter out entries where seqname is 'chr10'
filter_indices <- grepl("^chr", annotation_peak$seqname)
annotation_peak_filtered <- annotation_peak[filter_indices, ]

# Filter the rows in X
X_filtered <- X[filter_indices, ]

# Create the SummarizedExperiment object with the filtered data
atac <- SummarizedExperiment(assays = list(counts = X_filtered), 
                             rowRanges = GRanges(annotation_peak_filtered$seqname,
                             IRanges(annotation_peak_filtered$ranges)), 
                             colData = DataFrame(annotation_cells))
colnames(atac) <- annotation_cells$obs_id    

dim(atac) #peaks*cells

saveRDS(atac, paste0(par$temp_dir, "/scATAC/atac.rds"))



### Load RNA-seq and create sparsematrix
XX <- readMM(paste0(par$temp_dir, "/scRNA/X_matrix.mtx"))
XX <- t(XX)
annotation_gene <- read.csv(paste0(par$temp_dir, "/scRNA/annotation_gene.csv"), row.names = 1)
annotation_cells <- read.csv(paste0(par$temp_dir, "/scRNA/annotation_cells.csv"), row.names = 1)

rna <- as(XX, "CsparseMatrix")
rownames(rna) <- annotation_gene$location
colnames(rna) <- annotation_cells$obs_id

# Remove genes with zero expression across all cells
rna <- rna[Matrix::rowSums(rna)!=0,]

dim(rna) # genes*cells

saveRDS(rna, paste0(par$temp_dir, "/scRNA/rna.rds"))


## load cell topic probabilities and create cell-cluster matrix
cell_topic <- read.csv(paste0(par$cell_topic), row.names = 1)
print(dim(cell_topic))
# Derive cell kNN using this
cellkNN <- get.knn(cell_topic, k=n_topics)$nn.index
rownames(cellkNN) <- rownames(cell_topic)
print(dim(cellkNN))

saveRDS(cellkNN, paste0(par$temp_dir, "cellkNN.rds"))

