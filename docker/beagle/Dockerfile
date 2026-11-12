# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for running BEAGLE v5.5 (https://faculty.washington.edu/browning/beagle/beagle_5.5_17Dec24.pdf)

FROM ubuntu:24.04

LABEL maintainer="Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Install bcftools dependencies
RUN apt-get update && \
    apt-get install -y gcc make unzip wget bzip2 zlib1g-dev libbz2-dev liblzma-dev \
                       libcurl4-openssl-dev libssl-dev libncurses5-dev && \
    apt-get clean

# Install java
#ARG java=v8.511
RUN wget -O /tmp/jdk.tar.gz \
    https://download.oracle.com/java/27/latest/jdk-27_linux-x64_bin.tar.gz && \
    mkdir -p /opt/java && \
    tar -xzf /tmp/jdk.tar.gz -C /opt/java --strip-components=1 && \
    rm /tmp/jdk.tar.gz

ENV JAVA_HOME=/opt/java
ENV PATH="${JAVA_HOME}/bin:${PATH}"

# Download BEAGLE SOFTWARE and BREF3 software
RUN mkdir -p /opt/beagle && \
    wget -O /opt/beagle/beagle.jar \
    https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar && \
    wget -O /opt/beagle/bref3.jar \
    https://faculty.washington.edu/browning/beagle/bref3.27Feb25.75f.jar

# Launch bash at runtime
CMD ["/bin/bash"]
