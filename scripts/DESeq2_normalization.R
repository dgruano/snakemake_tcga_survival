library(DESeq2)
library(TCGAbiolinks)
library(ggplot2)
library(dplyr)

log_msg <- function(...) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] %s", ts, paste0(..., collapse = "")))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Expected 3 arguments: outfile, rds_file, biomart_file")
}
outfile <- args[1]
rds_file <- args[2]
biomart_file <- args[3]
outdir <- dirname(outfile)
GDC_dir <- dirname(outdir)

# create output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
log_msg("Output directory ready: ", outdir)
# set working directory to output directory
setwd(GDC_dir)
log_msg("Working directory set to: ", GDC_dir)

# load data
data = readRDS(rds_file)
log_msg("Loaded RDS file: ", rds_file)
log_msg("Initial dimensions (genes x samples): ", nrow(data), " x ", ncol(data))
# remove "Solid Tissue Normal" samples
samples_to_keep <- colData(data)$sample_type != "Solid Tissue Normal"
data <- data[, samples_to_keep]
log_msg("Samples kept after removing Solid Tissue Normal: ", sum(samples_to_keep),
                " / ", length(samples_to_keep))
log_msg("Dimensions after sample filtering (genes x samples): ", nrow(data), " x ", ncol(data))

if (ncol(data) == 0) {
    stop("No samples remain after filtering out Solid Tissue Normal")
}

# read biomaRt file
genes = read.csv(biomart_file)
log_msg("Loaded biomaRt file: ", biomart_file)
log_msg("biomaRt rows: ", nrow(genes))

required_cols <- c("ensembl_gene_id", "hgnc_symbol")
missing_cols <- setdiff(required_cols, colnames(genes))
if (length(missing_cols) > 0) {
    stop(paste0("biomaRt file is missing required columns: ",
                            paste(missing_cols, collapse = ", ")))
}

genes$ensembl_gene_id <- sub("\\..*$", "", genes$ensembl_gene_id)
# create dataframe with data rownames, and then remove everything after "."
data_rownames <- data.frame(
    ensembl_gene_id = sub("\\..*$", "", rownames(data)),
    row_idx = seq_len(nrow(data)),
    stringsAsFactors = FALSE
)
# merge with genes to get hgnc_symbol
data_merged <- merge(
    data_rownames,
    genes[, c("ensembl_gene_id", "hgnc_symbol")],
    by = "ensembl_gene_id"
)
log_msg("Rows after Ensembl merge: ", nrow(data_merged))

data_merged <- data_merged[!is.na(data_merged$hgnc_symbol) & data_merged$hgnc_symbol != "", ]
log_msg("Rows after removing empty HGNC symbols: ", nrow(data_merged))
# remove duplicates based on hgnc_symbol, keeping the first occurrence
data_merged <- data_merged[!duplicated(data_merged$hgnc_symbol), ]
log_msg("Rows after deduplicating HGNC symbols: ", nrow(data_merged))

# filter matrix by valid row indices from Ensembl ID merge
row_idx <- data_merged$row_idx
if (any(is.na(row_idx))) {
    log_msg("Warning: ", sum(is.na(row_idx)), " NA row indices found after merge; dropping them")
}
row_idx <- row_idx[!is.na(row_idx)]

if (length(row_idx) == 0) {
    stop("No genes left after Ensembl/HGNC mapping and filtering")
}

data_filtered <- data[row_idx, ]
# convert rownames to hgnc_symbol
rownames(data_filtered) <- data_merged$hgnc_symbol[!is.na(data_merged$row_idx)]
log_msg("Final dimensions for DESeq2 input (genes x samples): ",
                nrow(data_filtered), " x ", ncol(data_filtered))
# normalize using DESeq2
dds <- DESeqDataSetFromMatrix(countData = assay(data_filtered),
                              colData = colData(data),
                              design = ~1)
log_msg("DESeqDataSet created successfully")
# save normalized object
saveRDS(dds, file = outfile)
log_msg("Saved output RDS: ", outfile)