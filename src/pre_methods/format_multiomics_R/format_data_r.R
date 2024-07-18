library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(FigR)
library(optparse)

# Example usage
par <- list(
    temp_dir = 'output/figr',
    rna_rds = 'resources/grn-benchmark/multiomics_r/rna.rds',
    atac_rds = 'resources/grn-benchmark/multiomics_r/atac.rds'
)


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

saveRDS(atac, par$atac_rds)



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

saveRDS(rna, par$rna_rds)




