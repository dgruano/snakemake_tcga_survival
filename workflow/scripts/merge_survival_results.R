library(openxlsx)

# Setup logging helper
log_msg <- function(msg) {
  message("[", Sys.time(), "] ", msg)
}

log_msg("Starting merge_survival_results.R")

# Parse arguments
if (exists("snakemake")) {
  tcga_cohorts_str <- snakemake@params[["cohorts"]]
  output           <- snakemake@output[[1]]
} else {
  args             <- commandArgs(trailingOnly = TRUE)
  tcga_cohorts_str <- args[1]
  output           <- args[2]
}
if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, append = TRUE)
  sink(log_con, append = TRUE, type = "message")
}

log_msg(paste("TCGA cohorts (comma-separated):", tcga_cohorts_str))
log_msg(paste("Output directory:", output))

# Parse comma-separated cohort list
tryCatch({
  dir_list <- trimws(strsplit(tcga_cohorts_str, ",")[[1]])
  dir_list <- dir_list[dir_list != ""]
  log_msg(paste("Loaded", length(dir_list), "cohorts from comma-separated list"))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to parse TCGA cohorts:", conditionMessage(e)))
  stop(e)
})

out_path_survival <- dirname(output)
log_msg(paste("Working directory:", out_path_survival))

# set working directory
dir.create(out_path_survival, showWarnings = FALSE, recursive = TRUE)
setwd(out_path_survival)
log_msg("Working directory set and validated")

# helper function to merge survival tables
merge_survival_pvalues <- function(input_filename = "survival_pval.tsv",
                                   output_filename = "survival_pval_merged.xlsx") {

  log_msg(paste("Merging", input_filename, "-> ", output_filename))

  # create blank excel workbook
  OUT <- createWorkbook()
  sheets_added <- 0
  rows_total <- 0

  # loop through dirs to summarize the information
  for (dir in dir_list){
    file <- file.path(dir, input_filename)

    # Check if file exists
    if (!file.exists(file)) {
      log_msg(paste("WARNING: File not found for cohort", dir, "- skipping"))
      warning(paste("Missing file:", file))
      next
    }

    # Read data
    tryCatch({
      data <- read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
      log_msg(paste("  Read", nrow(data), "signatures from", dir))

      # add one sheet per cohort
      addWorksheet(OUT, dir)
      # write the survival data
      writeData(OUT, sheet = dir, x = data, rowNames = TRUE, colNames = TRUE)
      sheets_added <- sheets_added + 1
      rows_total <- rows_total + nrow(data)
    }, error = function(e) {
      log_msg(paste("ERROR: Failed to process cohort", dir, ":", conditionMessage(e)))
      warning(paste("Error processing", dir, ":", conditionMessage(e)))
    })
  }

  # save workbook
  tryCatch({
    if (file.exists(output_filename)) {
      log_msg(paste("WARNING: Overwriting existing file", output_filename))
      unlink(output_filename)
    }
    saveWorkbook(OUT, output_filename)
    log_msg(paste("Saved", output_filename, "with", sheets_added, "sheets,",
                  rows_total, "signatures"))
  }, error = function(e) {
    log_msg(paste("ERROR: Failed to save workbook", output_filename, ":", conditionMessage(e)))
    stop(e)
  })
}

log_msg("Starting merge operations")
overall_start <- Sys.time()

# merge survival pval tables
log_msg("Step 1: Merging survival p-values")
merge_survival_pvalues(input_filename = "survival_pval.tsv",
                       output_filename = "survival_pval_merged.xlsx")

# merge filtered tables
log_msg("Step 2: Merging filtered p-values")
merge_survival_pvalues(input_filename = "survival_pval_filtered.tsv",
                       output_filename = "survival_pval_filtered_merged.xlsx")

# save also a reorganized excel with the summary per gene/signature
log_msg("Step 3: Creating per-signature summary workbook")

# empty list to store each excel sheet data
data_list <- list()

# loop through tsvs
log_msg("  Collecting data from all cohorts...")
tryCatch({
  for (i in seq_along(dir_list)) {
    file_path <- file.path(dir_list[i], "survival_pval.tsv")
    # read tsv, skip if missing
    if (!file.exists(file_path)) {
      log_msg(paste("  WARNING: Missing file for", dir_list[i], "- skipping"))
      next
    }

    dat <- read.table(file_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  # save tsv in list
    data_list[[dir_list[i]]] <- dat
  }
  log_msg(paste("  Collected data from", length(data_list), "cohorts"))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to collect signature data:", conditionMessage(e)))
  stop(e)
})

# get gene/signature names
if (length(data_list) == 0) {
  log_msg("ERROR: No data collected from any cohort")
  stop("No valid data to merge")
}

all_signatures <- rownames(data_list[[1]])
log_msg(paste("  Found", length(all_signatures), "unique signatures"))

# create empty excel book
OUT <- createWorkbook()
signatures_processed <- 0

# build a table where rows = dirs and columns = percentiles for each signature
log_msg("  Building per-signature tables...")
for (signature in all_signatures) {
  # container for one feature across all directories
  df_row <- data.frame(matrix(nrow = length(data_list), ncol = ncol(data_list[[1]])))
  colnames(df_row) <- colnames(data_list[[1]])
  rownames(df_row) <- names(data_list)

  # fill from each dataset
  for (dir in names(data_list)) {
    df_row[dir, ] <- data_list[[dir]][signature, ]
  }

  # add worksheet
  addWorksheet(OUT, signature)
  writeData(OUT, sheet = signature, x = df_row,
            rowNames = TRUE, colNames = TRUE)
  signatures_processed <- signatures_processed + 1
}

# save workbook
tryCatch({
  if (file.exists("merged_per_signature.xlsx")) {
    log_msg("WARNING: Overwriting existing file merged_per_signature.xlsx")
    unlink("merged_per_signature.xlsx")
  }
  saveWorkbook(OUT, "merged_per_signature.xlsx")
  log_msg(paste("Saved merged_per_signature.xlsx with", signatures_processed,
                "signature sheets"))
}, error = function(e) {
  log_msg(paste("ERROR: Failed to save merged_per_signature.xlsx:", conditionMessage(e)))
  stop(e)
})

# Final summary
overall_elapsed <- difftime(Sys.time(), overall_start, units = "mins")
log_msg(paste("All merge operations completed successfully in", round(overall_elapsed, 2), "minutes"))
log_msg("merge_survival_results.R finished successfully")
