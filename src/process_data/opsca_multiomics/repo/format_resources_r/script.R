library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(FigR)

## VIASH START

par <- list(
    rna_matrix = 'output/scRNA/X_matrix.mtx',
    atac_matrix ='output/scATAC/X_matrix.mtx',
    rna_gene_annot = 'output/scRNA/annotation_gene.csv',
    rna_cell_annot = 'output/scRNA/annotation_cell.csv',
    atac_peak_annot = 'output/scATAC/annotation_gene.csv',
    atac_cell_annot = 'output/scATAC/annotation_cell.csv', 
    rna_rds = 'resources_test/grn-benchmark/multiomics_r/rna.rds',
    atac_rds = 'resources_test/grn-benchmark/multiomics_r/atac.rds'
)
## VIASH END


## Load atac-seq and create summarizedexperiment
X <- readMM(par$atac_matrix)
X <- t(X)
annotation_peak <- read.csv(par$atac_peak_annot, row.names = 1)
annotation_cells <- read.csv(par$atac_cell_annot, row.names = 1)

# Filter out entries where seqname is 'chr'
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
rownames(atac) <- paste(as.character(seqnames(atac)), as.character(ranges(atac)), sep=':')   

print(dim(atac)) #peaks*cells

saveRDS(atac, par$atac_rds)



### Load RNA-seq and create sparsematrix
XX <- readMM(par$rna_matrix)
XX <- t(XX)
annotation_gene <- read.csv(par$rna_gene_annot, row.names = 1)
annotation_cells <- read.csv(par$rna_cell_annot, row.names = 1)

rna <- as(XX, "CsparseMatrix")
rownames(rna) <- annotation_gene$location
colnames(rna) <- annotation_cells$obs_id

# Remove genes with zero expression across all cells
rna <- rna[Matrix::rowSums(rna)!=0,]

print(dim(rna)) # genes*cells

saveRDS(rna, par$rna_rds)