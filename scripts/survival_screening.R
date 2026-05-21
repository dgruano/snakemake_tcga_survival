library(DESeq2)
library(TCGAbiolinks)
library(ggplot2)
library(dplyr)
library(GSVA)
library(survminer)
library(survival)

log_msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste0(..., collapse = "")))
}

is_blank <- function(x) {
  is.na(x) | trimws(x) == ""
}

extract_genes <- function(x) {
  x <- x[!is_blank(x)]
  if (length(x) == 0) {
    return(character(0))
  }
  parts <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  unique(parts[parts != ""])
}

# define parameters
args <- commandArgs(trailingOnly = TRUE)
DPI <- args[1] %>% as.numeric()
THRESHOLD <- args[2] %>% as.numeric()
signatures_file <- args[3]
project <- args[4]
survival_table <- args[5]
GDC_dir <- survival_table %>% dirname() %>% dirname() %>% dirname() %>% dirname() 
# # correct aspect ratio of plots
# FACTOR = DPI / 100

# set working directory to output directory
setwd(GDC_dir)
# create PCA/project and survival/project directories
dir.create(paste0("./screening/PCA/", project), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0("./screening/survival/", project), showWarnings = FALSE,  recursive = FALSE)


#################################################################
###################### SURVIVAL SCREENING  ######################
#################################################################

# load normalized data
dds = readRDS(paste0("./DESeq2_normalized/", project,"_STAR_Counts_DESeq2.rds"))
# read gene_signatures.txt (tab separated)
gene_signatures <- read.delim(
  signatures_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fill = TRUE
)
log_msg("Loaded signatures file: ", signatures_file)
log_msg("Signature table dimensions (rows x cols): ", nrow(gene_signatures), " x ", ncol(gene_signatures))

if (ncol(gene_signatures) == 0) {
  stop("Signatures file has no columns")
}

# Some malformed files repeat the header as the first data row.
if (nrow(gene_signatures) > 0) {
  first_row <- as.character(gene_signatures[1, , drop = TRUE])
  if (length(first_row) == ncol(gene_signatures) &&
      all(trimws(first_row) == trimws(colnames(gene_signatures)))) {
    gene_signatures <- gene_signatures[-1, , drop = FALSE]
    log_msg("Detected duplicated header row in signatures file; first data row was removed")
  }
}

if (nrow(gene_signatures) == 0) {
  log_msg("No data rows in signatures file after cleanup; using column names as single-gene signatures")
}

signature_genes_raw <- lapply(gene_signatures, extract_genes)
if (nrow(gene_signatures) == 0) {
  signature_genes_raw <- as.list(colnames(gene_signatures))
  names(signature_genes_raw) <- colnames(gene_signatures)
}

signature_genes_raw <- signature_genes_raw[lengths(signature_genes_raw) > 0]
if (length(signature_genes_raw) == 0) {
  stop("No signatures with at least one gene were found in signatures file")
}

# extract normalized counts for ssGSEA and map signature IDs to expression row names
# in a case-insensitive way to avoid symbol case mismatches.
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
expr_genes <- rownames(norm_counts)
expr_gene_lookup <- setNames(expr_genes, toupper(expr_genes))
signature_genes <- lapply(signature_genes_raw, function(g) {
  mapped <- unname(expr_gene_lookup[toupper(g)])
  unique(mapped[!is.na(mapped)])
})

signature_sizes <- vapply(signature_genes, length, integer(1))
log_msg(
  "Parsed signatures: ",
  length(signature_genes),
  " total (",
  sum(signature_sizes == 0),
  " unmatched, ",
  sum(signature_sizes == 1),
  " single-gene, ",
  sum(signature_sizes >= 2),
  " multi-gene)"
)

unmatched_signatures <- names(signature_genes)[signature_sizes == 0]
if (length(unmatched_signatures) > 0) {
  preview <- paste(head(unmatched_signatures, 10), collapse = ", ")
  suffix <- ifelse(length(unmatched_signatures) > 10, " ...", "")
  log_msg(
    "Dropping signatures with no matching genes in expression matrix (n=",
    length(unmatched_signatures),
    "): ",
    preview,
    suffix
  )
}

signature_genes <- signature_genes[signature_sizes > 0]
if (length(signature_genes) == 0) {
  stop("No signatures could be matched to expression row names")
}

# load vector with percentiles
percentiles <- c(0.15, 0.20, 0.25, 0.33)

# classify signatures by number of matched genes
single_gene_signatures <- signature_genes[vapply(signature_genes, length, integer(1)) == 1]
multiple_gene_signatures <- signature_genes[vapply(signature_genes, length, integer(1)) >= 2]
# identify if any of the multiple gene signatures end with _UP or _DOWN and they have the same character string before that suffix
names_multiple <- names(multiple_gene_signatures)
base_names <- sub("(_UP|_DOWN)$", "", names_multiple)
duplicated_base_names <- base_names[duplicated(base_names)]

# perform ssGSEA only when >=1 valid multi-gene signature exists
if (length(multiple_gene_signatures) > 0) {
  scores <- gsva(
    as.matrix(norm_counts),
    multiple_gene_signatures,
    method = "ssgsea",
    ssgsea.norm = TRUE,
    verbose = FALSE
  )
  scores <- as.data.frame(t(scores))
  log_msg("Computed ssGSEA scores for multi-gene signatures: ", ncol(scores))
} else {
  scores <- data.frame(row.names = colnames(norm_counts))
  log_msg("No valid multi-gene signatures after ID matching; skipping ssGSEA step")
}

# combine _UP and _DOWN signatures with the same base name
for (base_name in duplicated_base_names) {
  up_name <- paste0(base_name, "_UP")
  down_name <- paste0(base_name, "_DOWN")
  if (up_name %in% colnames(scores) & down_name %in% colnames(scores)) {
    combined_scores <- scores[, up_name] - scores[, down_name]
    scores <- cbind(scores, combined_scores)
    colnames(scores)[ncol(scores)] <- paste0(base_name, "_COMBINED")
  }
}

# stratify using percentiles and save categorical variables in colData
for (signature in colnames(scores)) {
  for (percentile in percentiles) {
    threshold_low <- quantile(scores[, signature], probs = percentile)
    threshold_high <- quantile(scores[, signature], probs = 1 - percentile)
    categorical_var <- ifelse(scores[, signature] <= threshold_low, "Low",
                              ifelse(scores[, signature] >= threshold_high, "High", "Intermediate"))
    col_name <- paste0(signature, "_", percentile * 100, "pct")
    colData(dds)[, col_name] <- categorical_var
  }
}

# add the single_gene_signatures to the scores
for (signature in names(single_gene_signatures)) {
  gene_id <- single_gene_signatures[[signature]][1]
  scores[, signature] <- norm_counts[gene_id, ]
  for (percentile in percentiles) {
    threshold_low <- quantile(scores[, signature], probs = percentile)
    threshold_high <- quantile(scores[, signature], probs = 1 - percentile)
    categorical_var <- ifelse(scores[, signature] <= threshold_low, "Low",
                              ifelse(scores[, signature] >= threshold_high, "High", "Intermediate"))
    col_name <- paste0(signature, "_", percentile * 100, "pct")
    colData(dds)[, col_name] <- categorical_var
  }
}

# create empty p-value table
surv_pval_mat <- matrix(
  NA,
  nrow = length(colnames(scores)),
  ncol = length(percentiles),
  dimnames = list(colnames(scores), paste0(percentiles * 100, "pct"))
)

# loop through signatures and percentiles to plot PCA and survival
for (signature in colnames(scores)) {
  for (percentile in percentiles) {
    # prevent the loop from stopping if an error occurs
    tryCatch({
      col_name <- paste0(signature, "_", percentile * 100, "pct")
      dds_filtered <- colData(dds)
      # obtain the indices of non-Intermediate patients
      pat_to_remove <- which(dds_filtered[, col_name] != "Intermediate")
      # remove intermediate patients
      dds_filtered <- dds_filtered[pat_to_remove, ]
      # replace NA time to death with last follow-up time for censored cases
      notDead <- is.na(dds_filtered$days_to_death)
      dds_filtered$days_to_death[notDead] <- dds_filtered$days_to_last_follow_up[notDead]
      # create event column (s = TRUE if dead, FALSE if alive)
      dds_filtered$s <- grepl("dead|deceased", dds_filtered$vital_status, ignore.case = TRUE)
      # create grouping factor
      dds_filtered$type <- factor(dds_filtered[[col_name]], levels = c( "Low", "High"))
      # keep only required columns
      dds_filtered_surv <- dds_filtered[, c("days_to_death", "s", "type")]
      # convert days to months to ease readability
      dds_filtered_surv$days_to_death <- dds_filtered_surv$days_to_death / 365
      # survival model
      surv_obj <- Surv(time = as.numeric(dds_filtered_surv$days_to_death), event = dds_filtered_surv$s)
      fit <- survfit(surv_obj ~ type, data = dds_filtered_surv)
      # extract p-value
      # create the Kaplan-Meier plot only if p-value < THRESHOLD
      # create the PCA only under the same condition
      pval <- surv_pvalue(fit)$pval
      # save the p-value in the matrix
      surv_pval_mat[signature, paste0(percentile * 100, "pct")] <- pval
      if (pval < THRESHOLD) {
        # survival legend title
        if (pval < 0.00001) {
          surv_title <- paste0(project, " - ", signature, "\n", "Expression (p-val < 0.0001)")
        } else {
          surv_title <- paste0(project, " - ", signature, "\n", "Expression (p-val = ", round(pval, 5), ")")
        }
        # custom legend labels with sample sizes
        label.add.n <- function(x) {
          n <- sum(dds_filtered_surv$type == x)
          paste0(x, " (n=", n, ")")
        }
        legend_labels <- sapply(levels(dds_filtered_surv$type), label.add.n)

        # plot survival
        p <- ggsurvplot(
          fit,
          data = dds_filtered_surv,
          theme = theme_survminer(),
          risk.table = FALSE,
          pval = FALSE,
          conf.int = FALSE,
          fontsize = 50,
          xlab = "Time (years)",
          ylab = "Survival probability",
          legend.title = surv_title,
          legend.labs = legend_labels,
          palette = c("#255BA8", "#ED412B"),
          break.time.by = 1
        )
        ggsave(filename = paste0("./screening/survival/", project, "/Survival_", signature, "_", percentile * 100, "pct.png"), plot = p$plot, dpi = DPI, width = 6, height = 4)
        # filter dds object for PCA
        dds_filtered_pca <- dds[, rownames(dds_filtered)]
        # calculate PC1 and PC2 for PCA plot
        pca_data <- plotPCA(vst(dds_filtered_pca), intgroup = col_name, returnData = TRUE)
        percentVar <- round(100 * attr(pca_data, "percentVar"), 2)
        # plot PCA
        ggplot(pca_data, aes(PC1, PC2, color = get(col_name))) +
        geom_point(size = 3) +
        scale_color_manual(values = c("Low" = "#255BA8", "High" = "#ED412B")) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle(paste0(project, " - PCA - ", signature, " - ", percentile * 100, "pct")) +
        theme_minimal() +
        labs(color = "Expression groups")
        # save PCA as png
        ggsave(filename = paste0("./screening/PCA/", project, "/PCA_", signature, "_", percentile * 100, "pct.png"), dpi = DPI, width = 6, height = 4)
    }
    }, error = function(e) {

      # Print a warning but continue the loop
      message("⚠️ Error with signature '", signature,
              "', percentile ", percentile,
              ": ", conditionMessage(e))
      # store NA upon error
      surv_pval_mat[signature, paste0(percentile * 100, "pct")] <- NA
      # Continue to next iteration
      return(NULL)
    })
  }
}

# save scores matrix
write.table(
  scores,
  file = paste0("./screening/survival/", project,"/patient_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# save patient stratification on colData
tmp <- as.data.frame(colData(dds))
logical <- endsWith(colnames(tmp), "pct")
write.table(
  tmp[, logical],
  file = paste0("./screening/survival/", project,"/patient_scores_categorical.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# save survival p-value matrix
write.table(
  surv_pval_mat,
  file = paste0("./screening/survival/", project,"/survival_pval.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# create a filtered version of the p-value matrix
surv_pval_filtered <- surv_pval_mat
surv_pval_filtered[surv_pval_filtered >= THRESHOLD] <- ""

# save filtered p-value matrix
write.table(
  surv_pval_filtered,
  file = survival_table,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)