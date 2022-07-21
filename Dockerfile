FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y bedtools && \
    apt-get install -y build-essential wget \
    libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev \
    git tabix bc

## Code from: https://github.com/genome/docker-bcftools/blob/master/Dockerfile ##
RUN apt-get update && apt-get install -y \
  bzip2 \
  g++ \
  libbz2-dev \
  libcurl4-openssl-dev \
  liblzma-dev \
  make \
  ncurses-dev \
  wget \
  zlib1g-dev

ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV BCFTOOLS_VERSION=1.12

WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
  tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2

WORKDIR /tmp/bcftools-$BCFTOOLS_VERSION
RUN make prefix=$BCFTOOLS_INSTALL_DIR && \
  make prefix=$BCFTOOLS_INSTALL_DIR install

WORKDIR /
RUN ln -s $BCFTOOLS_INSTALL_DIR/bin/bcftools /usr/bin/bcftools && \
  rm -rf /tmp/bcftools-$BCFTOOLS_VERSION 

#####

RUN wget https://raw.githubusercontent.com/pj-sullivan/VCF_filter/main/VCF_filter.sh
RUN chmod +x VCF_filter.sh

RUN wget https://raw.githubusercontent.com/pj-sullivan/VCF_filter/main/gencode.v41.annotation.proteincoding.gtf.bed.gz
RUN wget https://raw.githubusercontent.com/pj-sullivan/VCF_filter/main/gencode.v41.annotation.proteincoding.omim.gtf.bed.gz
RUN wget https://raw.githubusercontent.com/pj-sullivan/VCF_filter/main/gencode.v41.annotation.omim.gtf.bed.gz
RUN wget https://raw.githubusercontent.com/pj-sullivan/VCF_filter/main/gencode.v41.annotation.gtf.bed.gz
