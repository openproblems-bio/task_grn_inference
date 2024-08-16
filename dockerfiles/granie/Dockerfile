# Use the base image
FROM  bioconductor/bioconductor_docker:devel-R-4.4.1


# Install required dependencies for the R packages
#RUN apt-get update && apt-get install -y \
#    libcurl4-openssl-dev \
#    libxml2-dev \
#    libssl-dev \
#    libcairo2-dev 
#    libxt-dev \
#    libopenblas-dev


#Install R packages
RUN R -e "install.packages(c('devtools'))"

# Set the working directory
WORKDIR /workspace

# Default command
CMD ["R"]

RUN R -e "remotes::install_version('Matrix', version = '1.6-3')"

#Version 1.9.4 of GRaNIE
RUN R -e "devtools::install_gitlab('grp-zaugg/GRaNIE@855bf3def33ad7af9353db79c7e21c9279035fb8', host = 'git.embl.de', subdir = 'src/GRaNIE', dependencies = TRUE, upgrade = 'never')"

#Version 0.2.1 of GRaNIEverse
RUN R -e "devtools::install_gitlab('grp-zaugg/GRaNIEverse@3224b042f25f0085eeaf9194042bcb89103ea962', host = 'git.embl.de', upgrade = 'never', dependencies = TRUE)"

#There seems to be an incompatibility with the Matrix and irlba package, downgrading Matrix seems to work

RUN R -e "install.packages('irlba',type='source')"

RUN R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86', 'EnsDb.Mmusculus.v79', 'BSgenome.Mmusculus.UCSC.mm39'))"



