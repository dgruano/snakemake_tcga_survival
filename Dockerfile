# syntax=docker/dockerfile:1
FROM rocker/r-ver:4.4.1

LABEL org.opencontainers.image.source="https://github.com/dgruano/snakemake_tcga_survival"
LABEL org.opencontainers.image.description="TCGA Survival Analysis Pipeline — R Bioconductor + Snakemake"
LABEL org.opencontainers.image.licenses="MIT"

# Install system deps for R packages, Snakemake, and bc
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    bc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Snakemake
RUN pip3 install --no-cache-dir snakemake==9.13.7

# Install R packages (BiocManager, then Bioc + CRAN packages)
RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'
RUN Rscript -e 'BiocManager::install(c( \
    "TCGAbiolinks", \
    "DESeq2", \
    "GSVA", \
    "biomaRt", \
    "survival", \
    "survminer", \
    "openxlsx" \
), update=FALSE, ask=FALSE, Ncpus=4)'

# Copy workflow files
COPY workflow/ /opt/workflow/workflow/
COPY config/ /opt/workflow/config/
COPY config.yaml /opt/workflow/config.yaml
COPY Snakefile /opt/workflow/Snakefile
COPY envs/ /opt/workflow/envs/
COPY profiles/ /opt/workflow/profiles/
COPY run_workflow.sh /usr/local/bin/run_tcga_survival

RUN chmod +x /usr/local/bin/run_tcga_survival

WORKDIR /opt/workflow
