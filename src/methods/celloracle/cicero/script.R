# library(Signac)
# library(Seurat)
# library(GenomicRanges)
# library(ggplot2)
# library(patchwork)
# library(SeuratWrappers)
library(cicero)
# library(AnnotationHub)

set.seed(2017)

data(cicero_data) # sample data. see the instruction how to build this from 10X data in the website
input_cds <- make_atac_cds(cicero_data, binarize = TRUE)

# - Preprocess the data
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

# *** if you are using Monocle 3 alpha, you need to run the following line as well!
# NOTE: Cicero does not yet support the Monocle 3 beta release (monocle3 package). We hope
# to update soon!
#input_cds <- preprocessCDS(input_cds, norm_method = "none")

input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                      reduction_method = 'tSNE', norm_method = "none")

print(input_cds)
tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

# - Run Cicero
data("human.hg19.genome")
sample_genome <- subset(human.hg19.genome, V1 == "chr18")
conns <- run_cicero(cicero_cds, sample_genome) # Takes a few minutes to run
print(head(conns))