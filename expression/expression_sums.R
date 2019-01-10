#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "dplyr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    message("Installing missing packages ...\n")
    tryCatch (silent = TRUE,
            install.packages(setdiff(packages, rownames(installed.packages())),
                             repos = "http://cran.us.r-project.org"),
            warning = function(bc) {
              source("http://bioconductor.org/biocLite.R")
              biocLite(setdiff(packages, rownames(installed.packages())))
            },
            error = function(bc) {
              source("http://bioconductor.org/biocLite.R")
              biocLite(setdiff(packages, rownames(installed.packages())))
            })
}

# Command parser
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
parser <- ArgumentParser(epilog = paste("Sums transcript-level TPM estimates",
                                        "to the gene level."))
parser$add_argument("input",
                    type = "character",
                    help = "input data file path")
parser$add_argument("biomart", 
                    type = "character", 
                    help = "biomaRt info file path")
parser$add_argument("output",
                    type = "character", 
                    help = "output data file path")
args <- parser$parse_args()

# Read input data
message("Reading data ...")
data <- read.table(args$input,
                   sep              = "\t", 
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Read biomaRt info
message("Reading biomaRt data ...")
ids <- read.table(args$biomart,
                  sep              = "\t",
                  header           = TRUE,
                  stringsAsFactors = FALSE)
ids <- ids[c("ensembl_transcript_id", "ensembl_gene_id")]
names(ids) <- c("ENSTID", "ENSGID")

# Add gene IDs to transcript data
message("Summing transcript TPM to gene level ...")
merged <- merge(data, ids, by = "ENSTID", all.x = TRUE)
merged$ENSTID <- NULL

# Remove NA-only rows
merged <- merged[!rowSums(is.na(merged)) > 0, ]

# Summarise
sums <- merged %>% group_by(ENSGID) %>% summarise_all(funs(sum(as.numeric(.))))

# Save to file
write.table(sums,
            file = args$output,
            sep = "\t",
            row.names = FALSE)
message("Done.")
