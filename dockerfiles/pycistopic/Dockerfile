# Use the base image
FROM janursa/scenicplus:19-08-2024

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
