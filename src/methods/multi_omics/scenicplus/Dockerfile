# Use the base image
FROM python:3.11.8

# Install scenicplus
RUN git clone https://github.com/aertslab/scenicplus
RUN pip install -e scenicplus[dev]
WORKDIR ..

# Install missing Python dependencies
# RUN pip install flatbuffers

# Install Tabix (from htslib)
RUN apt update && apt install -y tabix bedtools

# Install Java
RUN apt update
RUN apt install -y openjdk-17-jdk
ENV JAVA_HOME /usr/lib/jvm/java-17-openjdk-amd64/