suppressPackageStartupMessages({
  library(DESeq2)
  library(survival)
  library(survminer)
  library(ggplot2)
})

log_msg <- function(msg) message("[", Sys.time(), "] ", msg)

log_msg("===== survival_per_signature.R started =====")

if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, append = TRUE)
  sink(log_con, append = TRUE, type = "message")

  in_scores      <- snakemake@input[["scores"]]
  in_categorical <- snakemake@input[["categorical"]]
  in_clinical    <- snakemake@input[["clinical"]]
  in_deseq2      <- snakemake@input[["deseq2_file"]]
  out_pval       <- snakemake@output[[1]]
  project        <- snakemake@wildcards[["project"]]
  signature      <- snakemake@wildcards[["signature"]]
  THRESHOLD      <- snakemake@params[["THRESHOLD"]]
  DPI            <- snakemake@params[["DPI"]]
  percentiles_str <- snakemake@params[["PERCENTILES"]]
} else {
  args           <- commandArgs(trailingOnly = TRUE)
  in_scores      <- args[1]
  in_categorical <- args[2]
  in_clinical    <- args[3]
  in_deseq2      <- args[4]
  out_pval       <- args[5]
  project        <- args[6]
  signature      <- args[7]
  THRESHOLD      <- as.numeric(args[8])
  DPI            <- as.numeric(args[9])
  percentiles_str <- args[10]
}

THRESHOLD   <- as.numeric(THRESHOLD)
DPI         <- as.numeric(DPI)
percentiles <- as.numeric(strsplit(percentiles_str, ",")[[1]])

log_msg(paste("Project:", project))
log_msg(paste("Signature:", signature))
log_msg(paste("Percentiles:", paste(percentiles * 100, "pct", collapse = ", ")))

out_dir      <- dirname(out_pval)
survival_dir <- out_dir
pca_dir      <- file.path(dirname(dirname(out_dir)), "PCA", project)
dir.create(survival_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pca_dir,      showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load pre-computed inputs
# ---------------------------------------------------------------------------
scores      <- read.delim(in_scores,      row.names = 1, check.names = FALSE)
categorical <- read.delim(in_categorical, row.names = 1, check.names = FALSE)
clinical    <- read.delim(in_clinical,    row.names = 1, check.names = FALSE)
log_msg(paste("Loading DESeq2 object:", in_deseq2))
dds <- readRDS(in_deseq2)
log_msg(paste("Loaded:", nrow(dds), "genes x", ncol(dds), "samples"))

# ---------------------------------------------------------------------------
# Guard: signature may be absent if all genes were missing in this cohort
# ---------------------------------------------------------------------------
if (!(signature %in% colnames(scores))) {
  log_msg(paste("WARNING: signature", signature, "not found in scores; writing empty output"))
  empty <- matrix(
    NA_real_,
    nrow     = 1,
    ncol     = length(percentiles),
    dimnames = list(signature, paste0(percentiles * 100, "pct"))
  )
  write.table(empty, out_pval, sep = "\t", quote = FALSE, col.names = NA)
  log_msg("===== survival_per_signature.R completed (signature absent) =====")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Initialise p-value output matrix (1 row × P percentiles)
# ---------------------------------------------------------------------------
pval_mat <- matrix(
  NA_real_,
  nrow     = 1,
  ncol     = length(percentiles),
  dimnames = list(signature, paste0(percentiles * 100, "pct"))
)

plots_created <- 0

for (pct in percentiles) {
  col_name <- paste0(signature, "_", pct * 100, "pct")

  if (!(col_name %in% colnames(categorical))) {
    log_msg(paste("WARNING: column", col_name, "missing from categorical; skipping"))
    next
  }

  tryCatch({
    # Build survival dataframe for High / Low patients only
    grp       <- categorical[[col_name]]
    keep_idx  <- grp != "Intermediate"
    surv_df   <- data.frame(
      days      = as.numeric(clinical$days_to_death[keep_idx]),
      last_fu   = as.numeric(clinical$days_to_last_follow_up[keep_idx]),
      status    = clinical$vital_status[keep_idx],
      type      = factor(grp[keep_idx], levels = c("Low", "High")),
      row.names = rownames(clinical)[keep_idx]
    )

    # Censor: replace NA death time with last follow-up
    na_death           <- is.na(surv_df$days)
    surv_df$days[na_death] <- surv_df$last_fu[na_death]
    surv_df$event      <- grepl("dead|deceased", surv_df$status, ignore.case = TRUE)
    surv_df$time_years <- surv_df$days / 365

    surv_obj  <- Surv(time = surv_df$time_years, event = surv_df$event)
    fit       <- survfit(surv_obj ~ type, data = surv_df)
    cox_model <- coxph(surv_obj ~ type, data = surv_df)
    cox_sum   <- summary(cox_model)

    hr      <- round(cox_sum$coefficients[, "exp(coef)"],  2)
    hr_lo   <- round(cox_sum$conf.int[, "lower .95"],      2)
    hr_hi   <- round(cox_sum$conf.int[, "upper .95"],      2)
    pval    <- surv_pvalue(fit)$pval

    pval_mat[signature, paste0(pct * 100, "pct")] <- pval

    if (!is.na(pval) && pval < THRESHOLD) {
      plots_created <- plots_created + 1

      pval_text <- if (pval < 0.00001) "p < 0.00001" else paste0("p = ", signif(pval, 3))
      annotation_text <- paste0(pval_text, "\nHR = ", hr, " (", hr_lo, ", ", hr_hi, ")")
      surv_title <- paste0(
        project, "\n",
        signature, " (", pct * 100, "th percentile)\n",
        "HR = ", hr, " (95% CI ", hr_lo, "–", hr_hi, ")"
      )
      legend_labels <- sapply(levels(surv_df$type), function(x) {
        paste0(x, " (n=", sum(surv_df$type == x), ")")
      })

      # Kaplan-Meier plot
      p <- ggsurvplot(
        fit,
        data         = surv_df,
        theme        = theme_survminer(),
        risk.table   = FALSE,
        pval         = FALSE,
        conf.int     = FALSE,
        fontsize     = 5,
        xlab         = "Time (years)",
        ylab         = "Survival probability",
        title        = surv_title,
        legend.title = "Expression",
        legend.labs  = legend_labels,
        palette      = c("#255BA8", "#ED412B"),
        break.time.by = 3
      )
      p$plot <- p$plot +
        annotate("text", x = 0, y = 0.03,
                 label = annotation_text, hjust = 0, vjust = 0, size = 4)

      km_file <- file.path(survival_dir,
                           paste0("Survival_", signature, "_", pct * 100, "pct.png"))
      ggsave(km_file, plot = p$plot, dpi = DPI, width = 6, height = 5)
      log_msg(paste("  KM saved:", km_file, "| p =", signif(pval, 3), "| HR =", hr))

      # PCA — original approach: VST on the filtered patient subset
      dds_filtered_pca <- dds[, rownames(surv_df)]
      pca_data <- plotPCA(vst(dds_filtered_pca), intgroup = col_name, returnData = TRUE)
      pct_var  <- round(100 * attr(pca_data, "percentVar"), 2)

      ggplot(pca_data, aes(PC1, PC2, color = get(col_name))) +
        geom_point(size = 3) +
        scale_color_manual(values = c("Low" = "#255BA8", "High" = "#ED412B")) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        ggtitle(paste0(project, " - PCA - ", signature, " - ", pct * 100, "pct")) +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "white")) +
        labs(color = "Expression groups")

      pca_file <- file.path(pca_dir,
                            paste0("PCA_", signature, "_", pct * 100, "pct.png"))
      ggsave(pca_file, dpi = DPI, width = 6, height = 4)
      log_msg(paste("  PCA saved:", pca_file))
    }
  }, error = function(e) {
    log_msg(paste("ERROR at", pct * 100, "pct:", conditionMessage(e)))
    pval_mat[signature, paste0(pct * 100, "pct")] <<- NA_real_
  })
}

log_msg(paste("Plots created:", plots_created, "out of", length(percentiles), "percentiles"))

# Apply threshold filter (blank non-significant cells) before writing
pval_filtered        <- pval_mat
is_nonsig            <- is.na(pval_filtered) | pval_filtered >= THRESHOLD
pval_filtered[is_nonsig] <- ""

write.table(pval_filtered, out_pval, sep = "\t", quote = FALSE, col.names = NA)
log_msg(paste("Wrote", out_pval))

log_msg("===== survival_per_signature.R completed =====")
