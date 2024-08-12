
FROM nvidia/cuda:11.7.1-base-ubuntu20.04
# Install required dependencies for the R packages
RUN apt-get update && apt-get install -y \
    bedtools \
    build-essential \
    curl \
    bash \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install \
    scglue==0.3.2 \
    pyscenic==0.12.1 \
    numpy==1.23.4 \
    scanpy \
    networkx \
    pyarrow \
    cytoolz \
    scikit-misc \
    cuda-python

# Set the working directory
WORKDIR /workspace

# Default command
CMD ["python"]
