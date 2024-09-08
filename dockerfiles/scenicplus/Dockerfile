# Use the base image
FROM python:3.11.8

# Install scenicplus
RUN git clone https://github.com/aertslab/scenicplus
WORKDIR scenicplus/
RUN git checkout fa55dae55203951c0c42e359b65a28010cd544f6
RUN pip install .
WORKDIR ..

# Install HTSlib and tabix
RUN git clone https://github.com/samtools/htslib.git
WORKDIR htslib/
RUN autoreconf -i
RUN git submodule update --init --recursive
RUN ./configure
RUN make
RUN make install
WORKDIR ..

# Install Java
RUN apt update
RUN apt install -y openjdk-17-jdk
ENV JAVA_HOME /usr/lib/jvm/java-17-openjdk-amd64/
