library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)

## Step1: Peak-gene association testing

library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)

out_dir <- '../../output'

atac = readRDS(paste0(out_dir, "/scATAC/atac.rds"))
rna  = readRDS(paste0(out_dir, "/scRNA/rna.rds"))
cisCorr <- FigR::runGenePeakcorr(ATAC.se = atac,
                           RNAmat = rna,
                           genome = "hg38", # One of hg19, mm10 or hg38 
                           nCores = 40,
                           p.cut = NULL, # Set this to NULL and we can filter later
                           n_bg = 100)
write.csv(cisCorr, paste0(out_dir, "/infer/figr/grn/cisCorr.csv"), row.names = TRUE)

