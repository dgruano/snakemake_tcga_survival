library(biomaRt)
library(curl)

message("[", Sys.time(), "] Starting biomaRt gene annotation download")

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
outdir <- dirname(outfile)

message("[", Sys.time(), "] Output file: ", outfile)

# create output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
message("[", Sys.time(), "] Created output directory: ", outdir)

# Use the Ensembl Genes dataset
message("[", Sys.time(), "] Connecting to Ensembl database...")
mart <- NULL

tryCatch({
  mart <<- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  message("[", Sys.time(), "] Successfully connected to Ensembl")
}, error = function(e) {
  if (grepl("Multiple cache results", conditionMessage(e))) {
    message("[", Sys.time(), "] WARNING: biomaRt cache corrupted, clearing and retrying...")
    biomartCacheClear()
    tryCatch({
      mart <<- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
      message("[", Sys.time(), "] Successfully connected to Ensembl after cache clear")
    }, error = function(e2) {
      message("[", Sys.time(), "] ERROR: Failed to connect to Ensembl after cache clear")
      message("[", Sys.time(), "] Error message: ", conditionMessage(e2))
      stop(e2)
    })
  } else {
    message("[", Sys.time(), "] ERROR: Failed to connect to Ensembl")
    message("[", Sys.time(), "] Error message: ", conditionMessage(e))
    stop(e)
  }
})

# Retrieve mapping for protein-coding genes
message("[", Sys.time(), "] Retrieving protein-coding gene mappings from Ensembl...")
tryCatch({
  genes <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
    filters    = "biotype",
    values     = "protein_coding",
    mart       = mart
  )
  message("[", Sys.time(), "] Successfully retrieved ", nrow(genes), " protein-coding genes")
}, error = function(e) {
  message("[", Sys.time(), "] ERROR: Failed to retrieve gene mappings")
  message("[", Sys.time(), "] Error message: ", conditionMessage(e))
  stop(e)
})

# save file
message("[", Sys.time(), "] Writing results to file...")
tryCatch({
  write.csv(genes, file = outfile, row.names = FALSE)
  message("[", Sys.time(), "] Successfully wrote gene annotation file")
}, error = function(e) {
  message("[", Sys.time(), "] ERROR: Failed to write output file")
  message("[", Sys.time(), "] Error message: ", conditionMessage(e))
  stop(e)
})

message("[", Sys.time(), "] biomaRt download completed successfully")
