#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "clusterProfiler", "org.Hs.eg.db")
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

# Command parser --------------------------------------------------------------

suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser(epilog = paste("Performs GO and/or KEGG enrichment",
                                       "analyses of DEG data."))
parser$add_argument("input",
                    type    = "character",
                    help    = "input DEG file path")
parser$add_argument("biomart",
                    type    = "character",
                    help    = "path to biomaRt info file for gene IDs")
parser$add_argument("output",
                    type    = "character",
                    help    = "output file path")
parser$add_argument("-t", "--type",
                    type    = "character",
                    dest    = "type",
                    default = "GO",
                    metavar = "",
                    help    = "type [GO (default), GOslim, KEGG, KEGGM]")
args = parser$parse_args()

# Function definitions --------------------------------------------------------

# Function for reading biomaRt info
get_biomart_info <- function(biomart_file) {

    # Read info file
    message("Reading biomaRt info ...")
    info <- read.table(biomart_file,
                       sep              = "\t",
                       header           = TRUE,
                       stringsAsFactors = FALSE)

    # Remove unnecessary columns
    info <- info[c("ensemble_gene_id", "entrezgene")]

    # Convert to character
    info$entrezgene <- as.character(info$entrezgene)

    # Return biomaRt info
    return(info)
}

# Function for performing enrichment analyses
enrichment <- function(data, info, type = "GO") {
    
    # Perform specified enrichment
    if (type == "GO" | type == "GOslim") {

        # GO-enrichment
        enrich <- enrichGO(gene          = data$entrezgene,
                           universe      = info$entrezgene,
                           OrgDb         = org.Hs.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = args$p_cutoff,
                           qvalueCutoff  = args$q_cutoff,
                           readable      = TRUE)

        # GOslim-enrichment
        if (type == "GOslim") {

            # Simplify the GO-terms
            enrich <- simplify(enrich,
                               cutoff     = 0.7,
                               by         = "p.adjust",
                               select_fun = min)
        }
    } else if (type == "KEGG") {

        # KEGG-enrichment
        enrich <- enrichKEGG(gene         = data$entrezgene,
                             organism     = "hsa",
                             pvalueCutoff = args$p_cutoff)
    } else if (type == "KEGGM") {

        # KEGG Module-enrichment
        enrich <- enrichMKEGG(gene        = data$entrezgene,
                              organism    = "hsa")
    }

    # Return enrichment results
    return(enrich)
}

# Analysis --------------------------------------------------------------------

# Read DEG data
message("Reading DEG data ...")
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
data <- read.table(args$input,
                   sep              = "\t",
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Read biomaRt info
message("Reading biomaRt info ...")
info <- get_biomart_info(args$biomart)

# Perform enrichment analysis
message("Performing enrichment analysis ...")
enrich <- enrichment(data, info, args$type)

# Write results to file
write.table(enrich, args$output,
            sep       = "\t",
            row.names = FALSE)
message("Done.")
