FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ="America/New_York"

# Base OS deps (build tools + common libs) + libpng for matplotlib
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates git unzip zip \
    build-essential pkg-config \
    zlib1g-dev libbz2-dev liblzma-dev \
    libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev libsqlite3-dev \
    libfreetype6-dev libpng-dev \
    python3-venv python3-distutils python3-dev \
    # bio tools via apt instead of building
    samtools tabix \
    perl \
 && rm -rf /var/lib/apt/lists/*

# Install CPython 3.9.17 from source
ARG PYTHON_VERSION=3.9.17
RUN set -eux; \
    cd /tmp; \
    curl -fsSLO https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz; \
    tar -xzf Python-${PYTHON_VERSION}.tgz; \
    cd Python-${PYTHON_VERSION}; \
    ./configure --enable-optimizations; \
    make -j"$(nproc)"; \
    make install; \
    cd /; rm -rf /tmp/Python-${PYTHON_VERSION}*; \
    ln -s /usr/local/bin/python3 /usr/local/bin/python; \
    ln -s /usr/local/bin/pip3 /usr/local/bin/pip

# Make constraints global for all pip installs
COPY constraints.txt /tmp/constraints.txt
ENV PIP_CONSTRAINT=/tmp/constraints.txt \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_DEFAULT_TIMEOUT=180

# Clean any existing numpy/matplotlib remnants aggressively
RUN python - <<'PY'
import sys, site, pkgutil, shutil, pathlib
paths = set(site.getsitepackages() + [site.getusersitepackages()])
for p in list(paths):
    if not p: 
        continue
    for name in ("numpy", "matplotlib"):
        for m in pathlib.Path(p).glob(name):
            shutil.rmtree(m, ignore_errors=True)
        for m in pathlib.Path(p).glob(f"{name}-*.dist-info"):
            shutil.rmtree(m, ignore_errors=True)
        for m in pathlib.Path(p).glob(f"{name}-*.egg-info"):
            shutil.rmtree(m, ignore_errors=True)
print("Cleaned numpy/matplotlib from:", *paths, sep="\n  - ")
PY

# Install bedtools
RUN apt-get update && apt-get install -y --no-install-recommends bedtools \
    && rm -rf /var/lib/apt/lists/*

# Install tools + exact pins
RUN python -m pip install --no-cache-dir -U pip setuptools wheel \
 && pip install --no-cache-dir --upgrade --force-reinstall \
      "numpy==1.26.4" "matplotlib==3.8.4" "cython<3"

# Install MACS2. Build-from-source packages must reuse our pinned toolchain
RUN pip install --no-cache-dir --no-build-isolation MACS2==2.2.9.1

# Install Dictys without dependencies (we'll install them manually right after)
RUN pip install --no-cache-dir --no-build-isolation --no-deps \
  git+https://github.com/pinellolab/dictys.git@a82930fe8030af2785f9069ef5e909e49acc866f

# Install Dictys dependencies and more
RUN pip install --no-cache-dir --prefer-binary \
  pandas scipy networkx h5py threadpoolctl joblib \
  jupyter jupyterlab adjustText pyro-ppl docutils requests

# Install pyDNase and anndata without dependencies so it can't pin matplotlib<2
RUN pip install --no-cache-dir --no-build-isolation --no-deps pyDNase clint pysam packaging array_api_compat legacy-api-wrap zarr natsort anndata

# Install pybedtools version that works with cython<3
RUN pip install --no-cache-dir --no-build-isolation "pybedtools==0.9.1"

# Install pytorch
# RUN pip install --no-cache-dir --prefer-binary --index-url https://download.pytorch.org/whl/cpu torch

# HOMER prerequisites
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget perl unzip build-essential zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Install HOMER core + hg38 genome
RUN set -eux; \
    mkdir -p /opt/homer && cd /opt/homer; \
    curl -fsSLO http://homer.ucsd.edu/homer/configureHomer.pl; \
    chmod +x configureHomer.pl; \
    perl configureHomer.pl -install homer; \
    perl configureHomer.pl -install homerTools; \
    perl configureHomer.pl -install hg38
ENV PATH="/opt/homer/bin:${PATH}"

# hg38 annotations
RUN set -eux; \
    cd /opt/homer; \
    grep "hg38" update.txt > tmp.txt && mv tmp.txt update.txt; \
    cd update && ./updateUCSCGenomeAnnotations.pl ../update.txt

# Install CUDA
# RUN curl -fsSLo /etc/apt/preferences.d/cuda-repository-pin-600 \
#       https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin && \
#     curl -fsSLo /usr/share/keyrings/nvidia-cuda.gpg \
#       https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/3bf863cc.pub && \
#     echo "deb [signed-by=/usr/share/keyrings/nvidia-cuda.gpg] http://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/ /" \
#       > /etc/apt/sources.list.d/cuda.list && \
#     apt-get update && apt-get install -y --no-install-recommends cuda && \
#     rm -rf /var/lib/apt/lists/*

# Install AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install

CMD ["/bin/bash"]