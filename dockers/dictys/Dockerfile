FROM ubuntu:20.04

# Set non-interactive for apt
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV MPLCONFIGDIR=/tmp/matplotlib

# Install basic system packages
# Install basic system packages including wget
# Install wget and ca-certificates before downloading Miniforge
RUN apt update -y && apt install -y --no-install-recommends \
        wget \
        curl \
        ca-certificates \
        gcc \
        g++ \
        autoconf \
        automake \
        libtool \
        pkg-config \
        python3-dev \
        ninja-build \
        tzdata \
        libstdc++6 \
        cmake \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libcairo2-dev \
        libxt-dev \
        libopenblas-dev \
        bedtools && \
    rm -rf /var/lib/apt/lists/*

# Install Miniforge
RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3.sh -b -p /opt/conda && \
    rm Miniforge3.sh

# Create conda environment with packages
RUN /opt/conda/bin/conda init bash && \
    . /opt/conda/etc/profile.d/conda.sh && \
    /opt/conda/bin/mamba create -y -n dictys_env -c lingfeiwang -c conda-forge -c bioconda -c pytorch \
        python=3.10 \
        pip \
        dictys \
        pytorch \
        torchvision \
        torchaudio \
        cpuonly \
        jupyterlab \
        mudata \
		awscli \
		anndata \
        macs2 && \
    /opt/conda/bin/pip cache purge && \
    /opt/conda/bin/conda clean -a -y

# Default shell
SHELL ["/bin/bash", "-c"]

# Activate environment by default
SHELL ["/bin/bash", "-c"]
ENV PATH=/opt/conda/envs/dictys_env/bin:$PATH
CMD ["bash"]