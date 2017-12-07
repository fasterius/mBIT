#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages = c("argparse", "biomaRt")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  message("installing missing packages ...")
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
parser = ArgumentParser(epilog = "Download or update biomaRt info file.")
parser$add_argument("assembly",
                    type = "character",
                    help = "assembly to use [GRCh37 or GRCh38]")
parser$add_argument("output",
                    type = "character",
                    help = "output file path")
args <- parser$parse_args()

# Get assembly version
message("Getting assembly info ...")
suppressPackageStartupMessages(suppressWarnings(library("biomaRt")))
assembly = tolower(args$assembly)

# Get biomaRt information
if (assembly == "grch37") {

    # GRCh37 assembly
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  host    = "grch37.ensembl.org", 
                  path    = "/biomart/martservice",
                  dataset = "hsapiens_gene_ensembl")
} else if (assembly == "grch38") {

    # GRCh38 assembly
    mart <- useMart("ensembl",
                    dataset = "hsapiens_gene_ensembl")
} else {

    # Throw an error and stop
    stop(paste("incorrect assembly:", assembly))
}

# Get mart info and update info file
message("Getting gene/transcript info from biomaRt ...")
attributes <- c("ensembl_transcript_id",
                "ensembl_gene_id",
                "hgnc_symbol",
                "gene_biotype")
info <- getBM(attributes = attributes, mart = mart) 

# Save or update the stored assembly file
write.table(info, args$output, sep = "\t", quote = TRUE, row.names = FALSE)
