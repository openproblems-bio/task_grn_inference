# Use the base image
FROM ghcr.io/openproblems-bio/base_r:1.0.2

# Install required dependencies for the R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('dplyr', 'FNN', 'chromVAR', 'doParallel', 'BuenColors', 'FigR', 'BSgenome.Hsapiens.UCSC.hg38'), 
repos='http://cran.rstudio.com/')"

# Set the working directory
WORKDIR /workspace

# Default command
CMD ["R"]
# Use the base image
FROM ghcr.io/openproblems-bio/base_r:1.0.2

# Install required dependencies for the R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('dplyr', 'FNN', 'chromVAR', 'doParallel', 'BuenColors', 'FigR', 'BSgenome.Hsapiens.UCSC.hg38'), 
repos='http://cran.rstudio.com/')"

# Set the working directory
WORKDIR /workspace

# Default command
CMD ["R"]

