FROM lfwa/dictys-cpu:latest

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ="America/New_York"

# Install bio tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl unzip less groff \
    tabix samtools \
    && rm -rf /var/lib/apt/lists/*

# Install anndata
RUN pip install --no-cache-dir anndata

# Install AWS CLI v2
RUN curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip" \
    && unzip -q /tmp/awscliv2.zip -d /tmp \
    && /tmp/aws/install --update \
    && rm -rf /tmp/aws /tmp/awscliv2.zip

# Clean MACS2, install everything into dictys, add MACS2 wrapper
RUN source /opt/conda/etc/profile.d/conda.sh \
    && conda activate dictys \
    && (pip uninstall -y MACS2 || true) \
    && (conda remove -y macs2 || true) \
    && conda install -y -c conda-forge -c bioconda macs3 \
    && printf '#!/usr/bin/env bash\nexec macs3 "$@"\n' > "$(conda run -n dictys --no-capture-output python -c 'import sys,os; print(os.path.join(os.path.dirname(sys.executable), "macs2"))')" \
    && chmod +x "$(conda run -n dictys --no-capture-output python -c 'import sys,os; print(os.path.join(os.path.dirname(sys.executable), "macs2"))')" \
    && conda clean -afy

CMD ["/bin/bash"]