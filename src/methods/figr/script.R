library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)


## VIASH START
par <- list(
  multiomics_rna = "resources_test/grn-benchmark/multiomics_r/rna.rds",
  multiomics_atac = "resources_test/grn-benchmark/multiomics_r/atac.rds",
  temp_dir =  "output/temp_figr/",
  cell_topic = "resources/grn-benchmark/supp/cell_topic.csv",
  num_workers = 4,
  n_topics = 48,
  cisCorr = "output/figr/cisCorr.csv",
  prediction= "output/prediction.csv"
)
# meta <- list(
#   functionality_name = "my_method_r"
# )
## VIASH END

cellknn_func <- function(par) {
  ## load cell topic probabilities and create cell-cluster matrix
  cell_topic <- read.csv(paste0(par$cell_topic), row.names = 1)
  print(dim(cell_topic))
  # Derive cell kNN using this
  cellkNN <- get.knn(cell_topic, k=par$n_topics)$nn.index
  rownames(cellkNN) <- rownames(cell_topic)
  print(dim(cellkNN))
  saveRDS(cellkNN, paste0(par$temp_dir, "cellkNN.rds"))
}

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

## Step 2: create DORCs and smooth them 
dorc_genes_func <- functions(par){
  cisCorr = read.csv(paste0(par$temp_dir, "cisCorr.csv"))
  cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

  atac = readRDS(par$multiomics_atac)
  rna  = readRDS(par$multiomics_rna)

  allGenes = unique(cisCorr.filt$Gene) 
  dorcMat <- getDORCScores(ATAC.se = atac, # Has to be same SE as used in previous step
                          dorcTab = cisCorr.filt,
                          geneList = allGenes,
                          nCores = par$num_workers)

  cellkNN = readRDS(paste0(par$temp_dir, "cellkNN.rds"))
  # Smooth dorc scores using cell KNNs (k=n_topics)
  n_topics = par$n_topics
  dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:n_topics], mat = dorcMat, nCores = 4)

  # Smooth RNA using cell KNNs
  # This takes longer since it's all genes
  RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:n_topics], mat = rna,nCores = 4)

  # get peak gene connection
  peak_gene_figr = cisCorr.filt
  # peak_gene_figr_n = peak_gene_figr.groupby('Gene').apply(lambda df:df['PeakRanges'].shape[0])
  # np.max(peak_gene_figr_n.values), np.median(peak_gene_figr_n.values)
  # print('In the peak-gene associations: number of  CIS ', peak_gene_figr.PeakRanges.unique().shape[0], ', gene ', peak_gene_figr.Gene.unique().shape[0])
  # print('number of DORC genes ', (peak_gene_figr_n.values >= 10).sum())
  # peak_gene_figr = peak_gene_figr[['PeakRanges', 'Gene']]
  # peak_gene_figr.columns = ['peak','target']
  # peak_gene_figr.to_csv(f'{out_dir}/infer/figr/peak_gene.csv')

  write.csv(cisCorr.filt, paste0(par$temp_dir, "cisCorr.filt.csv"))
  saveRDS(RNAmat.s, paste0(par$temp_dir, "RNAmat.s.RDS"))
  saveRDS(dorcMat.s, paste0(par$temp_dir, "dorcMat.s.RDS"))
}
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

  # write.csv(figR.d, paste0(out_dir, "/infer/figr/grn/figR.d.csv"))
}



# filter based on enrichment 
figr_grn = figr_grn[figr_grn['Enrichment.P']<0.05]
# filter bsaed on correlatoon
figr_grn = figr_grn[figr_grn['Corr.P']<0.05]
# filter thoes that are 0 score 
figr_grn = figr_grn[figr_grn.Score!=0]
# subset columns
figr_grn = figr_grn[['Motif', 'DORC', 'Score']]
figr_grn = figr_grn.reset_index(drop=True)
figr_grn.columns = ['source', 'target','weight']
figr_grn.to_csv(f'{out_dir}/infer/figr/grn/figr_grn.csv')