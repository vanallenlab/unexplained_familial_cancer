# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for vanallenlab/bcftools

FROM us.gcr.io/google.com/cloudsdktool/google-cloud-cli:588.0.0

LABEL maintainer="Noah Fields <Noah_Fields@dfci.harvard.edu>"

# Install bcftools dependencies
RUN apt-get update && \
    apt-get install -y gcc make unzip wget bzip2 zlib1g-dev libbz2-dev liblzma-dev \
                       libcurl4-openssl-dev libssl-dev libncurses5-dev && \
    apt-get clean

# Install bcftools
#ARG BCFTOOLS_VERISON=1.19
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xvf bcftools-1.19.tar.bz2 && \
    cd bcftools-1.19 && \
    ./configure --prefix=/ && \
    make && \
    make install

# Install PLINK 1.9
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210614.zip && \
    unzip plink_linux_x86_64_20210614.zip -d /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink && \
    rm plink_linux_x86_64_20210614.zip

RUN rm -r bcftools-1.19
RUN rm bcftools-1.19.tar.bz2

# Launch bash at runtime
CMD ["/bin/bash"]
