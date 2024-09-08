# Use the base image
FROM  bioconductor/bioconductor_docker:devel-R-4.4.1


#Install R packages
RUN R -e "install.packages(c('devtools'))"

# Set the working directory
WORKDIR /workspace

# Default command
CMD ["R"]

#Version 1.9.4 of GRaNIE.
RUN R -e "devtools::install_gitlab('grp-zaugg/GRaNIE@v1.9.4', host = 'git.embl.de', subdir = 'src/GRaNIE', dependencies = TRUE, upgrade = 'never')"

#Version 0.3.1 of GRaNIEverse
RUN R -e "devtools::install_gitlab('grp-zaugg/GRaNIEverse@v1.3.1', host = 'git.embl.de', upgrade = 'never', dependencies = TRUE)"

#There is an incompatibility with the Matrix and irlba package, see https://github.com/satijalab/seurat/issues/8000
RUN R -e "install.packages('irlba',type='source', rebuild = TRUE)"

# biovizbase is only needed to create the ATAC assay along with annotation
RUN R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86', 'EnsDb.Mmusculus.v79', 'BSgenome.Mmusculus.UCSC.mm39', 'biovizBase', 'Signac', 'glmGamPoi'))"
