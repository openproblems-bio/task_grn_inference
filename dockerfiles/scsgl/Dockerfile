# Use the base image
FROM python:3.8

# Install dependencies
RUN apt-get update
RUN apt-get install -y r-base time
RUN Rscript -e "install.packages('pcaPP')"

# Install scSGL
RUN git clone https://github.com/Single-Cell-Graph-Learning/scSGL
WORKDIR "/scSGL"
ENV SCSGL_PATH $PWD
