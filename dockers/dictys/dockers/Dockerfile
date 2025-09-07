# Start with Miniconda base image
FROM continuumio/miniconda3:23.10.0-1

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Copy environment file
COPY environment.yml environment.yml

# Create conda environment
RUN conda env create -f environment.yml && \
    conda clean -afy

# Set environment path
SHELL ["conda", "run", "-n", "dictys", "/bin/bash", "-c"]
ENV PATH /opt/conda/envs/dictys/bin:$PATH

# (Optional) Set default env activation
RUN echo "conda activate dictys" >> ~/.bashrc

# Default command
CMD ["python"]