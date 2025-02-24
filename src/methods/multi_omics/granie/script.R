set.seed(42)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Matrix))
library(GRaNIEverse)
library(GRaNIE)
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm39))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SummarizedExperiment))


## VIASH START
par <- list(
  multiomics_rna_r = "resources_test/grn_benchmark/inference_data/op_rna.rds",
  multiomics_atac_r = "resources_test/grn-benchmark/multiomics_atac.rds",
  preprocessing_clusteringMethod = 1, # Seurat::FindClusters: (1 = original Louvain algorithm, 2 = Louvain algorithm with multilevel refinement, 3 = SLM algorithm, 4 = Leiden algorithm)
  preprocessing_clusterResolution = 14, # Typically between 5 and 20
  preprocessing_RNA_nDimensions = 50, # Default 50
  genomeAssembly = "hg38",
  GRaNIE_corMethod = "spearman",
  GRaNIE_includeSexChr = TRUE,
  GRaNIE_promoterRange = 250000,
  GRaNIE_TF_peak_fdr_threshold = 0.2,
  GGRaNIE_peak_gene_fdr_threshold = 0.2,
  num_workers = 4,
  peak_gene = "output/granie/peak_gene.csv", # not yet implemented, should I?
  prediction= "output/granie/prediction.csv",
  useWeightingLinks = FALSE,
  forceRerun = FALSE
)
## VIASH END


cat("Content of par list:")
str(par)

#### STANDARD ASSIGNMENTS ###
file_seurat = "seurat_granie.qs"
outputDir = dirname(par$prediction)

if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}


#########################
# Downloading resources #
#########################
file_hocomoco_v12 = "https://s3.embl.de/zaugg-web/GRaNIE/TFBS/hg38/PWMScan_HOCOMOCOv12_H12INVIVO.tar.gz"
destfile <- "PWMScan_HOCOMOCOv12_H12INVIVO.tar.gz"
if (!file.exists(destfile)) {
  options(timeout = 1200)
  download.file(file_hocomoco_v12, destfile)
}
# Define the directory to extract the files to
exdir <- "PWMScan_HOCOMOCOv12_H12INVIVO"
GRaNIE_TFBSFolder = paste0(exdir, "/H12INVIVO")
if (!file.exists(GRaNIE_TFBSFolder)) {
  untar(destfile, exdir = exdir)
}

if (par$genomeAssembly == "hg38"){
  file_RNA_URL = "https://s3.embl.de/zaugg-web/GRaNIEverse/features_RNA_hg38.tsv.gz"

} else if (par$genomeAssembly == "mm10") {
  file_RNA_URL = "https://s3.embl.de/zaugg-web/GRaNIEverse/features_RNA_mm10.tsv.gz"
}

file_RNA <- paste0("features_RNA_", par$genomeAssembly, ".tsv.gz")
if (!file.exists(file_RNA)) {
  options(timeout = 1200)
  download.file(file_RNA_URL, file_RNA)
}

print('Donwnloading resources finished')

###################
# Preprocess data #
###################

if (par$forceRerun | !file.exists(file_seurat)) {
 # Sparse matrix
 rna.m = readRDS(par$multiomics_rna_r)
 seurat_object <- CreateSeuratObject(count = rna.m, project = "PBMC", min.cells = 1, min.features = 1, assay = "RNA")
 
 # RangedSummarizedExperiment
 atac = readRDS(par$multiomics_atac_r)
 
 # Extract counts and metadata from the RangedSummarizedExperiment
  atac_counts <- assays(atac)$counts
  
  rownames(atac_counts) =  paste0(seqnames(rowRanges(atac)) %>% as.character(), ":", start(rowRanges(atac)), "-", end(rowRanges(atac)))
  
  # Create a ChromatinAssay
  chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = NULL,
   min.cells = 1,
   min.features = 1,
   colData = DataFrame(colData(atac))
  )
 
  if (par$genomeAssembly == "hg38"){
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    
  } else if (par$genomeAssembly == "mm10") {
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  }
  
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- par$genomeAssembly
  Annotation(chrom_assay) <- annotations
 
  # Unify cells
  # Identify the common cells between the RNA and ATAC assays
  common_cells <- intersect(colnames(seurat_object[["RNA"]]), colnames(chrom_assay))
  
  # Subset the Seurat object to include only the common cells
  chrom_assay <- subset(chrom_assay, cells = common_cells)

  seurat_object[["peaks"]] = chrom_assay
    
  qs::qsave(seurat_object, "seurat_granie.qs")
    
} else {

  seurat_object = qs::qread(file_seurat)
  
}

output_seuratProcessed = paste0(outputDir, "/seuratObject.qs")

###################
# Preprocess data #
###################

# Take output from preprocessing steps
GRaNIE_file_peaks = paste0(outputDir, "/atac.pseudobulkFromClusters_res", par$preprocessing_clusterResolution, "_mean.tsv.gz")
GRaNIE_file_rna = paste0(outputDir, "/rna.pseudobulkFromClusters_res", par$preprocessing_clusterResolution, "_mean.tsv.gz")
GRaNIE_file_metadata = paste0(outputDir, "/metadata_res", par$preprocessing_clusterResolution, "_mean.tsv.gz")

if (file.exists(GRaNIE_file_peaks) & file.exists(GRaNIE_file_metadata) & file.exists(GRaNIE_file_rna) & !par$forceRerun) {
  
  cat("Preprocessing skipped because all files already exist anf forceRerun = FALSE.")
  
} else {
  seurat_object = prepareSeuratData_GRaNIE(seurat_object, 
                                           outputDir = outputDir,
                                           saveSeuratObject = TRUE,
                                           genome = par$genomeAssembly,
                                           assayName_RNA = "RNA", normRNA = "SCT", nDimensions_RNA = par$preprocessing_RNA_nDimensions, recalculateVariableFeatures = NULL,
                                           assayName_ATAC_raw = "peaks", 
                                           normATAC = "LSI", LSI_featureCutoff = "q0", nDimensions_ATAC = 50, dimensionsToIgnore_LSI_ATAC = 1,
                                           integrationMethod = "WNN", WNN_knn = 20,
                                           pseudobulk_source = "cluster",
                                           countAggregation = "mean",
                                           clusteringAlgorithm = par$preprocessing_clusteringMethod, 
                                           clusterResolutions = par$preprocessing_clusterResolution,
                                           minCellsPerCluster = 25,
                                           forceRerun = FALSE
      )
  
}



##############
# Run GRaNIE #
##############

GRN = runGRaNIE(
  dir_output = outputDir,
  datasetName = "undescribed",
  GRaNIE_file_peaks,
  GRaNIE_file_rna,
  GRaNIE_file_metadata,
  TFBS_source = "custom",
  TFBS_folder = GRaNIE_TFBSFolder,
  genomeAssembly = par$genomeAssembly,
  normalization_peaks = "none",
  idColumn_peaks = "peakID",
  normalization_rna = "none",
  idColumn_RNA = "ENSEMBL",
  includeSexChr = par$GRaNIE_includeSexChr,
  minCV = 0,
  minNormalizedMean_peaks = NULL,
  minNormalizedMean_RNA = NULL,
  minSizePeaks = 5,
  corMethod = par$GRaNIE_corMethod,
  promoterRange = par$GRaNIE_promoterRange,
  useGCCorrection = FALSE,
  TF_peak.fdr.threshold = par$GRaNIE_TF_peak_fdr_threshold,
  peak_gene.fdr.threshold = par$GRaNIE_peak_gene_fdr_threshold,
  runTFClassification = FALSE,
  runNetworkAnalyses = FALSE,
  nCores = par$num_workers,
  forceRerun = TRUE
)

# Post-process GRN
connections.df = getGRNConnections(GRN, 
                                   include_TF_gene_correlations = TRUE, 
                                   include_peakMetadata = TRUE, 
                                   include_TFMetadata = TRUE, 
                                   include_geneMetadata = TRUE)

final.df = connections.df %>%
  dplyr::select(TF.name, gene.name, TF_gene.r) %>%
  dplyr::rename(source = TF.name, target = gene.name)

if (par$useWeightingLinks) {
  final.df = dplyr::mutate(final.df, weight = abs(.data$TF_gene.r))
} else {
  final.df = dplyr::mutate(final.df, weight = 1)
}

final.df %>%
  dplyr::select(source, target, weight) %>%
  readr::write_csv(par$prediction)

