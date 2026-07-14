# syntax=docker/dockerfile:1
FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.source="https://github.com/dgruano/snakemake_tcga_survival"
LABEL org.opencontainers.image.description="TCGA Survival Analysis Pipeline — R Bioconductor + Snakemake"
LABEL org.opencontainers.image.licenses="MIT"

# Install Snakemake and bc (for percentile math)
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-pip \
    python3-dev \
    bc \
    && pip3 install --no-cache-dir snakemake==9.13.7 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install additional R packages not in the base image
RUN Rscript -e ' \
    BiocManager::install(c( \
        "TCGAbiolinks", \
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
