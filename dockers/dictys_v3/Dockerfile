FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ="America/New_York"

RUN apt update && \
    apt install -y \
    wget build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev \
    libreadline-dev libffi-dev libsqlite3-dev libbz2-dev liblzma-dev \
    git zip wget software-properties-common \
    # For Dictys
    libfreetype-dev libfreetype6 libfreetype6-dev pkg-config && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

ARG PYTHON_VERSION=3.9.17

RUN cd /tmp && \
    wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz && \
    tar -xvf Python-${PYTHON_VERSION}.tgz && \
    cd Python-${PYTHON_VERSION} && \
    ./configure --enable-optimizations && \
    make && make install && \
    cd .. && rm Python-${PYTHON_VERSION}.tgz && rm -r Python-${PYTHON_VERSION} && \
    ln -s /usr/local/bin/python3 /usr/local/bin/python && \
    ln -s /usr/local/bin/pip3 /usr/local/bin/pip && \
    python -m pip install --upgrade pip && \
    rm -r /root/.cache/pip
					
RUN pip install --upgrade pip && \
    pip install --no-cache-dir \
    jupyterlab MACS2
    

# Samtools
RUN cd /usr/bin && \
    wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar -vxjf samtools-1.17.tar.bz2 && \
    cd samtools-1.17 && \
    ./configure --prefix=/usr/bin/samtools && \
    make && \
    make install

# Downgrade numpy for Dictys
RUN pip install numpy==1.23.5

# Dictys
RUN pip install git+https://github.com/pinellolab/dictys.git


# Install HOMER
RUN mkdir /opt/homer && \
    cd /opt/homer && \
    wget http://homer.ucsd.edu/homer/configureHomer.pl && \
    perl configureHomer.pl -install homer

RUN cd /opt/homer && \
    cat update.txt | grep "hg38" > tmp.txt && \
    mv tmp.txt update.txt && \
    cd update && \
    ./updateUCSCGenomeAnnotations.pl ../update.txt

RUN wget -O /etc/apt/preferences.d/cuda-repository-pin-600 https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin && \
    apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/3bf863cc.pub && \
    add-apt-repository "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/ /"

RUN apt update && \
    apt install -y cuda
    
ENV PATH="${PATH}:/opt/homer/bin"
ENV PATH="/usr/bin/samtools/bin:${PATH}"

CMD ["/bin/bash"]