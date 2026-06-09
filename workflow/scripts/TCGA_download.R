suppressPackageStartupMessages({
    library(TCGAbiolinks)
})

if (exists("snakemake")) {
  project <- snakemake@wildcards[["cohort"]]
  outfile <- snakemake@output[[1]]
} else {
  args    <- commandArgs(trailingOnly = TRUE)
  project <- args[1]
  outfile <- args[2]
}
if (exists("snakemake")) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con, append = TRUE)
  sink(log_con, append = TRUE, type = "message")
}
outdir <- dirname(outfile)
root_dir <- dirname(outdir)

# create rds output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# set working directory to root of output directory
setwd(root_dir)

# create query for project
query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

# download the project with retry logic
GDCdownload(query = query,
            method = "api",
            files.per.chunk = NULL)
data <- GDCprepare(query = query)

# save the data as .rds file
cat("Saving data for project", project, "to", outfile, "\n")
saveRDS(data, file = outfile)


# clean up GDC cache
gdc_temp_dir <- file.path(root_dir, "GDCdata", project)
if (dir.exists(gdc_temp_dir)) {
    cat("Cleaning up GDC cache at", gdc_temp_dir, "\n")
    unlink(gdc_temp_dir, recursive = TRUE, force = TRUE)
}
