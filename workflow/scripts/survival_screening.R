suppressPackageStartupMessages({
  library(DESeq2)
  library(TCGAbiolinks)
  library(ggplot2)
  library(dplyr)
  library(GSVA)
  library(survminer)
  library(survival)
})


# Logging helper function
log_msg <- function(msg) {
  message("[", Sys.time(), "] ", msg)
}

log_msg("===== survival_screening.R started =====")

# define parameters
if (exists("snakemake")) {
  DPI             <- snakemake@params[["DPI"]]
  THRESHOLD       <- snakemake@params[["THRESHOLD"]]
  signatures_file <- snakemake@input[["signatures_file"]]
  project         <- snakemake@wildcards[["project"]]
  survival_table  <- snakemake@output[[1]]
  percentiles_str <- snakemake@params[["PERCENTILES"]]
} else {
  args            <- commandArgs(trailingOnly = TRUE)
  DPI             <- args[1]
  THRESHOLD       <- args[2]
  signatures_file <- args[3]
  project         <- args[4]
  survival_table  <- args[5]
  percentiles_str <- args[6]
}
if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, append = TRUE)
  sink(log_con, append = TRUE, type = "message")
}
DPI <- as.numeric(DPI)
THRESHOLD <- as.numeric(THRESHOLD)
GDC_dir <- survival_table %>% dirname() %>% dirname() %>% dirname() %>% dirname()

# Parse percentiles from comma-separated string
percentiles <- as.numeric(strsplit(percentiles_str, ",")[[1]])

# Log command-line arguments
log_msg(paste("DPI:", DPI))
log_msg(paste("THRESHOLD:", THRESHOLD))
log_msg(paste("Signatures file:", signatures_file))
log_msg(paste("Project:", project))
log_msg(paste("Survival table output:", survival_table))
log_msg(paste("Percentiles:", paste(percentiles * 100, "pct", collapse = ", ")))
log_msg(paste("GDC directory:", GDC_dir))

# Validate input files exist
if (!file.exists(signatures_file)) {
  log_msg(paste("ERROR: Signatures file not found:", signatures_file))
  stop("Signatures file not found: ", signatures_file)
}
log_msg("Input files validated")

# set working directory to output directory
setwd(GDC_dir)
log_msg(paste("Working directory set to:", getwd()))

# create PCA/project and survival/project directories
pca_dir <- paste0("./screening/PCA/", project)
survival_dir <- paste0("./screening/survival/", project)
dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(survival_dir, showWarnings = FALSE, recursive = FALSE)
log_msg(paste("Output directories created:", pca_dir, "and", survival_dir))


#################################################################
###################### SURVIVAL SCREENING  ######################
#################################################################

log_msg("Starting survival screening pipeline...")
start_time <- Sys.time()

# load normalized data
log_msg(paste("Loading DESeq2 object from:", paste0("./DESeq2_normalized/", project,"_STAR_Counts_DESeq2.rds")))
tryCatch({
  dds <- readRDS(paste0("./DESeq2_normalized/", project,"_STAR_Counts_DESeq2.rds"))
  log_msg(paste("Loaded DESeq2 object with", nrow(dds), "genes x", ncol(dds), "samples"))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to load DESeq2 object:", conditionMessage(e)))
  stop(e)
})

# read gene_signatures.txt (tab separated)
log_msg(paste("Loading gene signatures from:", signatures_file))
tryCatch({
  gene_signatures <- read.delim(signatures_file, header = TRUE, stringsAsFactors = FALSE)
  log_msg(paste("Loaded", length(gene_signatures), "signature columns"))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to load gene signatures:", conditionMessage(e)))
  stop(e)
})

# Percentiles already loaded from command-line argument

# classify signatures in < 2 genes and >= 2 genes (column values)
log_msg("Classifying signatures...")
single_gene_signatures <- gene_signatures[sapply(gene_signatures, function(x) length(unlist(strsplit(x, ","))) < 2)]
multiple_gene_signatures <- gene_signatures[sapply(gene_signatures, function(x) length(unlist(strsplit(x, ","))) >= 2)]
log_msg(paste("Found", length(single_gene_signatures), "single-gene signatures and", length(multiple_gene_signatures), "multi-gene signatures"))

# identify if any of the multiple gene signatures end with _UP or _DOWN and they have the same character string before that suffix
names_multiple <- names(multiple_gene_signatures)
base_names <- sub("(_UP|_DOWN)$", "", names_multiple)
duplicated_base_names <- base_names[duplicated(base_names)]
if (length(duplicated_base_names) > 0) {
  log_msg(paste("Found", length(duplicated_base_names), "paired UP/DOWN signatures to combine:", paste(duplicated_base_names, collapse = ", ")))
}

# extract normlalized counts for ssGSEA
log_msg("Estimating size factors and extracting normalized counts...")
tryCatch({
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  log_msg(paste("Extracted normalized counts:", nrow(norm_counts), "genes x", ncol(norm_counts), "samples"))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to estimate size factors:", conditionMessage(e)))
  stop(e)
})

# perform ssgsea for multiple gene signatures
if (length(multiple_gene_signatures) > 0 && ncol(as.data.frame(multiple_gene_signatures)) > 0) {
  log_msg("Performing ssGSEA for multi-gene signatures...")
  tryCatch({
    scores <- gsva(norm_counts, as.list(multiple_gene_signatures), method = "ssgsea", ssgsea.norm = TRUE)
    scores <- as.data.frame(t(scores))
    log_msg(paste("ssGSEA completed: generated scores for", ncol(scores), "signatures"))
  }, error = function(e) {
    log_msg(paste("ERROR: ssGSEA failed:", conditionMessage(e)))
    stop(e)
  })
} else {
  log_msg("No multi-gene signatures found; creating empty scores dataframe")
  scores <- data.frame(row.names = colnames(norm_counts))
}

# combine _UP and _DOWN signatures with the same base name
if (length(duplicated_base_names) > 0) {
  log_msg("Combining paired UP/DOWN signatures...")
  for (base_name in duplicated_base_names) {
    up_name <- paste0(base_name, "_UP")
    down_name <- paste0(base_name, "_DOWN")
    if (up_name %in% colnames(scores) & down_name %in% colnames(scores)) {
      combined_scores <- scores[, up_name] - scores[, down_name]
      scores <- cbind(scores, combined_scores)
      colnames(scores)[ncol(scores)] <- paste0(base_name, "_COMBINED")
      log_msg(paste("  Created combined score:", paste0(base_name, "_COMBINED")))
    }
  }
}

# stratify using percentiles and save categorical variables in colData
log_msg("Stratifying multi-gene signatures by percentiles...")
stratification_count <- 0
for (signature in colnames(scores)) {
  for (percentile in percentiles) {
    threshold_low <- quantile(scores[, signature], probs = percentile)
    threshold_high <- quantile(scores[, signature], probs = 1 - percentile)
    categorical_var <- ifelse(scores[, signature] <= threshold_low, "Low",
                              ifelse(scores[, signature] >= threshold_high, "High", "Intermediate"))
    col_name <- paste0(signature, "_", percentile * 100, "pct")
    colData(dds)[, col_name] <- categorical_var
    stratification_count <- stratification_count + 1
  }
}
log_msg(paste("Created", stratification_count, "stratified variables for multi-gene signatures"))

# add the single_gene_signatures to the scores
log_msg(paste("Processing", length(single_gene_signatures), "single-gene signatures..."))
for (signature in names(single_gene_signatures)) {
  tryCatch({
    # check if the gene is present in the normalized counts, skip if not
    if (!(signature %in% rownames(norm_counts))) {
      log_msg(paste("WARNING: Gene", signature, "not found in normalized counts; skipping"))
      next
    }
    scores[, signature] <- norm_counts[signature, ]
    for (percentile in percentiles) {
      threshold_low <- quantile(scores[, signature], probs = percentile)
      threshold_high <- quantile(scores[, signature], probs = 1 - percentile)
      categorical_var <- ifelse(scores[, signature] <= threshold_low, "Low",
                                ifelse(scores[, signature] >= threshold_high, "High", "Intermediate"))
      col_name <- paste0(signature, "_", percentile * 100, "pct")
      colData(dds)[, col_name] <- categorical_var
      stratification_count <- stratification_count + 1
    }
  }, error = function(e) {
    log_msg(paste("WARNING: Failed to process single-gene signature", signature, ":", conditionMessage(e)))
  })
}
log_msg(paste("Total stratified variables created:", stratification_count))

# create empty p-value table
surv_pval_mat <- matrix(
  NA,
  nrow = length(colnames(scores)),
  ncol = length(percentiles),
  dimnames = list(colnames(scores), paste0(percentiles * 100, "pct"))
)
log_msg(paste("Initialized p-value matrix:", nrow(surv_pval_mat), "signatures x", ncol(surv_pval_mat), "percentiles"))

# loop through signatures and percentiles to plot PCA and survival
log_msg("Starting survival analysis and plot generation...")
total_combos <- length(colnames(scores)) * length(percentiles)
current_combo <- 0
plots_created <- 0

for (signature in colnames(scores)) {
  for (percentile in percentiles) {
    current_combo <- current_combo + 1
    # prevent the loop from stopping if an error occurs
    tryCatch({
      col_name <- paste0(signature, "_", percentile * 100, "pct")
      dds_filtered <- colData(dds)
      # obtain the indices of non-Intermediate patients
      pat_to_remove <- which(dds_filtered[, col_name] != "Intermediate")
      # remove intermediate patients
      dds_filtered <- dds_filtered[pat_to_remove, ]
      low_count <- sum(dds_filtered[[col_name]] == "Low")
      high_count <- sum(dds_filtered[[col_name]] == "High")

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
      # Cox proportional hazards model for HR
      cox_model <- coxph(surv_obj ~ type, data = dds_filtered_surv)
      cox_summary <- summary(cox_model)

      # Extract HR (High vs Low)
      hr <- round(cox_summary$coefficients[,"exp(coef)"], 2)
      hr_conf_low <- round(cox_summary$conf.int[,"lower .95"], 2)
      hr_conf_high <- round(cox_summary$conf.int[,"upper .95"], 2)

      # extract p-value
      # create the Kaplan-Meier plot only if p-value < THRESHOLD
      # create the PCA only under the same condition
      pval <- surv_pvalue(fit)$pval
      # save the p-value in the matrix
      surv_pval_mat[signature, paste0(percentile * 100, "pct")] <- pval

      if (pval < THRESHOLD) {
        plots_created <- plots_created + 1
        # format p-value text
        if (pval < 0.00001) {
          pval_text <- "p < 0.00001"
        } else {
          pval_text <- paste0("p = ", signif(pval, 3))
        }

        # format HR and p-value for the figure bottom left information
        annotation_text <- paste0(
          pval_text, "\n",
          "HR = ", hr, " (", hr_conf_low, ", ", hr_conf_high, ")"
        )

        # format figure title
        surv_title <- paste0(
          project, "\n",
          signature, " (", percentile * 100, "th percentile)\n",
          "HR = ", hr, " (95% CI ", hr_conf_low, "–", hr_conf_high, ")"
        )

        # custom legend labels with sample sizes
        label.add.n <- function(x) {
          n <- sum(dds_filtered_surv$type == x)
          paste0(x, " (n=", n, ")")
        }
        legend_labels <- sapply(levels(dds_filtered_surv$type), label.add.n)

        # create survival plot
        p <- ggsurvplot(
          fit,
          data = dds_filtered_surv,
          theme = theme_survminer(),
          risk.table = FALSE,
          pval = FALSE,
          conf.int = FALSE,
          fontsize = 5,
          xlab = "Time (years)",
          ylab = "Survival probability",
          title = surv_title,
          legend.title = "Expression",
          legend.labs = legend_labels,
          palette = c("#255BA8", "#ED412B"),
          break.time.by = 3
        )

        # include p-value in survival figure
        p$plot <- p$plot +
          annotate(
            "text",
            x = 0,
            y = 0.03,
            label = annotation_text,
            hjust = 0,
            vjust = 0,
            size = 4
          )

        # save file
        survival_plot_file <- paste0("./screening/survival/", project, "/Survival_", signature, "_", percentile * 100, "pct.png")
        ggsave(filename = survival_plot_file, plot = p$plot, dpi = DPI, width = 6, height = 5)
        log_msg(paste("  Created survival plot for", signature, "at", percentile * 100, "pct (Low: n=", low_count, ", High: n=", high_count, "); p=", signif(pval, 3), "; HR=", hr, ")"))

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
        theme(plot.background = element_rect(fill = "white")) +
        labs(color = "Expression groups")
        # save PCA as png
        pca_plot_file <- paste0("./screening/PCA/", project, "/PCA_", signature, "_", percentile * 100, "pct.png")
        ggsave(filename = pca_plot_file, dpi = DPI, width = 6, height = 4)
    }
    }, error = function(e) {
      log_msg(paste("ERROR: Analysis failed for signature", signature, "at", percentile * 100, "pct:", conditionMessage(e)))
      # store NA upon error
      surv_pval_mat[signature, paste0(percentile * 100, "pct")] <- NA
      # Continue to next iteration
      return(NULL)
    })
  }
}

log_msg(paste("Analysis complete: created", plots_created, "significant plots out of", total_combos, "analyses"))

# save scores matrix
log_msg("Writing output files...")
scores_file <- paste0("./screening/survival/", project, "/patient_scores.tsv")
tryCatch({
  write.table(
    scores,
    file = scores_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  log_msg(paste("Wrote patient scores:", nrow(scores), "samples x", ncol(scores), "signatures to", scores_file))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to write patient scores:", conditionMessage(e)))
  stop(e)
})

# save patient stratification on colData
categorical_file <- paste0("./screening/survival/", project, "/patient_scores_categorical.tsv")
tryCatch({
  tmp <- as.data.frame(colData(dds))
  logical <- endsWith(colnames(tmp), "pct")
  write.table(
    tmp[, logical],
    file = categorical_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  log_msg(paste("Wrote categorical scores:", nrow(tmp), "samples x", sum(logical), "stratified variables to", categorical_file))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to write categorical scores:", conditionMessage(e)))
  stop(e)
})

# save survival p-value matrix
pval_file <- paste0("./screening/survival/", project, "/survival_pval.tsv")
tryCatch({
  write.table(
    surv_pval_mat,
    file = pval_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  log_msg(paste("Wrote p-value matrix to", pval_file))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to write p-value matrix:", conditionMessage(e)))
  stop(e)
})

# create a filtered version of the p-value matrix
surv_pval_filtered <- surv_pval_mat
surv_pval_filtered[surv_pval_filtered >= THRESHOLD] <- ""

# save filtered p-value matrix
tryCatch({
  write.table(
    surv_pval_filtered,
    file = survival_table,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  log_msg(paste("Wrote filtered p-value matrix (threshold < ", THRESHOLD, ") to", survival_table))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to write filtered p-value matrix:", conditionMessage(e)))
  stop(e)
})

# Summary and completion
elapsed <- difftime(Sys.time(), start_time, units = "mins")
log_msg(paste("===== survival_screening.R completed successfully in", round(elapsed, 2), "minutes =====" ))
