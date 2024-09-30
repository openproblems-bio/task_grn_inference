# Use the base image
FROM python:3.11.8

# Install scenicplus
RUN git clone https://github.com/aertslab/scenicplus
WORKDIR scenicplus/
RUN git checkout fa55dae55203951c0c42e359b65a28010cd544f6
RUN pip install .
WORKDIR ..

# Install missing Python dependencies
RUN pip install flatbuffers

# Install Tabix (from htslib)
RUN apt update && apt install -y tabix bedtools

# Install Java
RUN apt update
RUN apt install -y openjdk-17-jdk
ENV JAVA_HOME /usr/lib/jvm/java-17-openjdk-amd64/
