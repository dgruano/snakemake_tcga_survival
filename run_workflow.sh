#!/bin/bash
# run_workflow.sh — wrapper for the TCGA Survival Snakemake pipeline
# Called from the Galaxy UDT. Galaxy's CWD is the working directory.
set -euo pipefail

GALAXY_CWD="$(pwd)"
TCGA_PROJECTS_FILE="$1"
SIGNATURES_FILE="$2"
PERCENTILES_IN="$3"
THRESHOLD="$4"
DPI="$5"
CORES="${GALAXY_SLOTS:-4}"

WORKFLOW_DIR="/opt/workflow"
cd "$WORKFLOW_DIR"

echo "=== TCGA Survival Pipeline ==="
echo "Galaxy CWD:     $GALAXY_CWD"
echo "Projects file:  $TCGA_PROJECTS_FILE"
echo "Signatures file: $SIGNATURES_FILE"
echo "Percentiles:    $PERCENTILES_IN"
echo "Threshold:      $THRESHOLD"
echo "DPI:            $DPI"
echo "Cores:          $CORES"

# Place user inputs where Snakemake expects them
cp "$TCGA_PROJECTS_FILE" config/tcga_projects.txt
cp "$SIGNATURES_FILE" config/gene_signatures.txt

# Convert comma-separated percentiles to Snakemake list-of-floats format
# e.g. "5,10,25,50" -> "[0.05,0.10,0.25,0.50]"
PCT_LIST="["
IFS=',' read -ra PCTS <<< "$PERCENTILES_IN"
FIRST=true
for p in "${PCTS[@]}"; do
    p="$(echo "$p" | xargs)"  # trim whitespace
    DEC=$(echo "scale=2; $p/100" | bc -l)
    if [ "$FIRST" = true ]; then
        PCT_LIST+="$DEC"
        FIRST=false
    else
        PCT_LIST+=",$DEC"
    fi
done
PCT_LIST+="]"
echo "Percentiles as decimals: $PCT_LIST"

# Create output directories
mkdir -p results logs

# Run Snakemake (no-use-conda — all R/Bioc packages pre-installed in container)
snakemake --cores "$CORES" \
    --config DPI="$DPI" \
             THRESHOLD="$THRESHOLD" \
             PERCENTILES="$PCT_LIST" \
             pathvars:results="$WORKFLOW_DIR/results" \
             pathvars:logs="$WORKFLOW_DIR/logs" \
    --latency-wait 60 \
    --no-use-conda \
    -p 2>&1

echo "=== Workflow complete. Zipping results... ==="

# Zip results back into the Galaxy working directory
cd "$WORKFLOW_DIR/results"
zip -r "$GALAXY_CWD/tcga_survival_results.zip" . 2>&1

echo "=== Done ==="
ls -lh "$GALAXY_CWD/tcga_survival_results.zip"
