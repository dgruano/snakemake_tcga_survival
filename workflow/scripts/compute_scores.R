suppressPackageStartupMessages({
  library(DESeq2)
  library(GSVA)
  library(BiocParallel)
})

log_msg <- function(msg) message("[", Sys.time(), "] ", msg)

log_msg("===== compute_scores.R started =====")

if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, append = TRUE)
  sink(log_con, append = TRUE, type = "message")

  deseq2_file     <- snakemake@input[["deseq2_file"]]
  signatures_file <- snakemake@input[["signatures_file"]]
  out_scores      <- snakemake@output[["scores"]]
  out_categorical <- snakemake@output[["categorical"]]
  out_clinical    <- snakemake@output[["clinical"]]
  percentiles_str <- snakemake@params[["PERCENTILES"]]
  n_threads       <- snakemake@threads
} else {
  args            <- commandArgs(trailingOnly = TRUE)
  deseq2_file     <- args[1]
  signatures_file <- args[2]
  out_scores      <- args[3]
  out_categorical <- args[4]
  out_clinical    <- args[5]
  percentiles_str <- args[6]
  n_threads       <- 1L
}

percentiles <- as.numeric(strsplit(percentiles_str, ",")[[1]])
log_msg(paste("Percentiles:", paste(percentiles * 100, "pct", collapse = ", ")))
log_msg(paste("Threads:", n_threads))

dir.create(dirname(out_scores), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load DESeq2 object (size factors already estimated by DESeq2_normalization.R)
# ---------------------------------------------------------------------------
log_msg(paste("Loading DESeq2 object:", deseq2_file))
dds <- readRDS(deseq2_file)
log_msg(paste("Loaded:", nrow(dds), "genes x", ncol(dds), "samples"))

norm_counts <- counts(dds, normalized = TRUE)
log_msg("Extracted normalised counts")

# ---------------------------------------------------------------------------
# Read and classify signatures
# ---------------------------------------------------------------------------
gene_signatures <- read.delim(signatures_file, header = TRUE, stringsAsFactors = FALSE)
log_msg(paste("Loaded", ncol(gene_signatures), "signature columns"))

is_single <- sapply(gene_signatures, function(col) {
  length(na.omit(col[col != ""])) < 2
})
single_gene_signatures   <- gene_signatures[,  is_single, drop = FALSE]
multiple_gene_signatures <- gene_signatures[, !is_single, drop = FALSE]
log_msg(paste(
  "Single-gene:", ncol(single_gene_signatures),
  "| Multi-gene:", ncol(multiple_gene_signatures)
))

# ---------------------------------------------------------------------------
# ssGSEA for multi-gene signatures
# ---------------------------------------------------------------------------
if (ncol(multiple_gene_signatures) > 0) {
  log_msg(paste("Running ssGSEA with", n_threads, "thread(s)..."))
  t0 <- proc.time()
  scores <- gsva(
    norm_counts,
    as.list(multiple_gene_signatures),
    method    = "ssgsea",
    ssgsea.norm = TRUE,
    BPPARAM   = MulticoreParam(n_threads)
  )
  scores <- as.data.frame(t(scores))
  log_msg(paste(
    "ssGSEA done in", round((proc.time() - t0)["elapsed"], 1), "s;",
    ncol(scores), "signatures"
  ))
} else {
  log_msg("No multi-gene signatures; creating empty scores dataframe")
  scores <- data.frame(row.names = colnames(norm_counts))
}

# ---------------------------------------------------------------------------
# Combine paired _UP / _DOWN signatures → _COMBINED
# ---------------------------------------------------------------------------
names_multi   <- colnames(scores)
base_names    <- sub("(_UP|_DOWN)$", "", names_multi)
paired_bases  <- unique(base_names[duplicated(base_names)])

for (base in paired_bases) {
  up   <- paste0(base, "_UP")
  down <- paste0(base, "_DOWN")
  if (up %in% colnames(scores) && down %in% colnames(scores)) {
    scores[[paste0(base, "_COMBINED")]] <- scores[[up]] - scores[[down]]
    log_msg(paste("Created combined score:", paste0(base, "_COMBINED")))
  }
}

# ---------------------------------------------------------------------------
# Append single-gene scores (raw normalised expression)
# ---------------------------------------------------------------------------
for (sig in colnames(single_gene_signatures)) {
  if (!(sig %in% rownames(norm_counts))) {
    log_msg(paste("WARNING: gene", sig, "not found in normalised counts; skipping"))
    next
  }
  scores[[sig]] <- norm_counts[sig, rownames(scores)]
}
log_msg(paste("Total signatures in score matrix:", ncol(scores)))

# ---------------------------------------------------------------------------
# Stratify every signature × percentile → High / Intermediate / Low
# ---------------------------------------------------------------------------
categorical <- data.frame(row.names = rownames(scores))

for (sig in colnames(scores)) {
  for (pct in percentiles) {
    lo  <- quantile(scores[[sig]], probs = pct,     na.rm = TRUE)
    hi  <- quantile(scores[[sig]], probs = 1 - pct, na.rm = TRUE)
    grp <- ifelse(scores[[sig]] <= lo, "Low",
           ifelse(scores[[sig]] >= hi, "High", "Intermediate"))
    categorical[[paste0(sig, "_", pct * 100, "pct")]] <- grp
  }
}
log_msg(paste("Stratification complete:", ncol(categorical), "columns"))

# ---------------------------------------------------------------------------
# Export clinical columns needed by survival_per_signature.R
# ---------------------------------------------------------------------------
clinical_cols <- c("days_to_death", "days_to_last_follow_up", "vital_status")
clinical      <- as.data.frame(colData(dds))[, clinical_cols, drop = FALSE]
log_msg("Extracted clinical data")

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
write.table(scores,      out_scores,      sep = "\t", quote = FALSE, col.names = NA)
log_msg(paste("Wrote patient_scores.tsv:", nrow(scores), "x", ncol(scores)))

write.table(categorical, out_categorical, sep = "\t", quote = FALSE, col.names = NA)
log_msg(paste("Wrote patient_scores_categorical.tsv:", nrow(categorical), "x", ncol(categorical)))

write.table(clinical,    out_clinical,    sep = "\t", quote = FALSE, col.names = NA)
log_msg(paste("Wrote clinical_data.tsv:", nrow(clinical), "x", ncol(clinical)))

log_msg("===== compute_scores.R completed =====")
