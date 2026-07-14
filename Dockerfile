# syntax=docker/dockerfile:1
FROM rocker/r-ver:4.4.1

LABEL org.opencontainers.image.source="https://github.com/dgruano/snakemake_tcga_survival"
LABEL org.opencontainers.image.description="TCGA Survival Analysis Pipeline — R Bioconductor + Snakemake"
LABEL org.opencontainers.image.licenses="MIT"

# Install system deps for R packages + bc
RUN apt-get update && apt-get install -y --no-install-recommends \
    bc \
    curl \
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

# Install Miniforge (conda) for Snakemake (v9.x from conda-forge)
RUN curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    -o /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH="/opt/conda/bin:${PATH}"

# Install Snakemake via conda-forge
RUN mamba install -c conda-forge -y snakemake && \
    mamba clean -afy

# Configure R to use Posit binary package manager (avoids compiling from source)
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))' \
    >> /usr/local/lib/R/etc/Rprofile.site

# Install R packages (binary from Posit PPM + BiocManager for Bioc packages)
RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'
RUN Rscript -e 'Sys.setenv(BIOCONDUCTOR_USE_CONTAINER_REPOSITORY="TRUE"); \
    BiocManager::install(c( \
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
