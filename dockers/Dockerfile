# syntax=docker/dockerfile:1

FROM continuumio/miniconda3:latest

# set up working directory
WORKDIR /workspace

# install mamba for faster solves
RUN conda install -y -c conda-forge mamba

# create environment with base deps
RUN mamba create -y -n dictys \
    -c lingfeiwang -c bioconda -c conda-forge -c pytorch -c nvidia \
    python=3.10 dictys pytorch torchvision torchaudio pytorch-cuda=11.7 mudata pip

# activate environment and install pip-only deps (torch/cu118 wheels)
RUN /bin/bash -c "source activate dictys && \
    pip install torch torchvision torchaudio \
    --index-url https://download.pytorch.org/whl/cu118"

# make environment default
SHELL ["conda", "run", "-n", "dictys", "/bin/bash", "-c"]

# optional: verify installation
RUN python -c "import torch; print('Torch CUDA available:', torch.cuda.is_available())"

ENTRYPOINT ["/bin/bash"]
