#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "clusterProfiler", "org.Hs.eg.db", "ggplot2")
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
                                       "analyses of gene or DEG data."))
parser$add_argument("input",
                    type    = "character",
                    help    = "input DEG file path")
parser$add_argument("biomart",
                    type    = "character",
                    help    = "path to biomaRt info file for gene IDs")
parser$add_argument("output",
                    type    = "character",
                    help    = "output figure file path")
parser$add_argument("-o", "--output-list",
                    action  = "store_true",
                    dest    = "output_list",
                    help    = "also output the full term list")
parser$add_argument("-t", "--type",
                    type    = "character",
                    dest    = "type",
                    default = "GO",
                    metavar = "",
                    help    = "type: GO [default], GOslim, KEGG, KEGGM")
parser$add_argument("-f", "--FDR-cutoff",
                    type    = "double",
                    dest    = "fdr_cutoff",
                    default = 0.01,
                    metavar = "",
                    help    = "FDR cutoff [default: 0.01]")
parser$add_argument("-n", "--top-n-terms",
                    type    = "integer",
                    dest    = "top_n_terms",
                    default = 10,
                    metavar = "",
                    help    = "top n terms to plot [default: 10]")
args = parser$parse_args()

# Function definitions --------------------------------------------------------

# Function for reading biomaRt info
get_biomart_info <- function(biomart_file) {

    # Read info file
    info <- read.table(biomart_file,
                       sep              = "\t",
                       header           = TRUE,
                       stringsAsFactors = FALSE)

    # Remove unnecessary columns, duplicates and missing values
    info <- unique(info[c("ensembl_gene_id", "entrezgene")])
    info <- info[!is.na(info$entrezgene), ]

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
                           keyType       = "ENTREZID",
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = args$fdr_cutoff,
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
                             pvalueCutoff = args$fdr_cutoff)
    } else if (type == "KEGGM") {

        # KEGG Module-enrichment
        enrich <- enrichMKEGG(gene        = data$entrezgene,
                              organism    = "hsa")
    }

    # Return enrichment results
    return(enrich)
}


# Function for plotting the top enriched terms
plot_top_terms <- function(enrich, output, top_n_terms = 10) {

    # Sort and get top terms
    enrich <- enrich[order(enrich$p.adjust), ]
    data <- head(enrich, top_n_terms)

    # Separate long terms with newlines
    data$Description <- gsub("([^ ]+ [^ ]+) ", "\\1\n", data$Description)

    # Set factors for plotting order
    data$Description <- factor(data$Description,
                               levels = rev(data$Description))

    # Plot
    gg <- ggplot(data, aes(x    = Description,
                           y    = -log10(p.adjust),
                           fill = -log10(p.adjust))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_classic() +
        labs(title = paste("Top", top_n_terms, "enriched", args$type, "terms"),
             x     = NULL,
             y     = expression(log[10](FDR))) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_fill_gradient(low   = "#a6c6f2",
                            high  = "#0d2d59",
                            guide = FALSE)

    # Save to output
    ggsave(output, gg, dpi = 300, height = 7, width = 7)
}

# Analysis --------------------------------------------------------------------

# Read DEG data
message("Reading gene data ...")
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("ggplot2"))
data <- read.table(args$input,
                   sep              = "\t",
                   header           = TRUE,
                   quote            = "",
                   stringsAsFactors = FALSE)

# Read biomaRt info
message("Reading biomaRt info ...")
info <- get_biomart_info(args$biomart)

# Merge data with biomaRt info to get Entrez IDs
message("Merging with gene data ...")
data <- merge(data, info, by.x = "ENSGID", by.y = "ensembl_gene_id")

# Perform enrichment analysis
message("Performing enrichment analysis ...")
enrich <- enrichment(data, info, args$type)
enrich <- head(enrich, nrow(enrich))

# Plot and save to file
plot_top_terms(enrich,
               args$output,
               args$top_n_terms)

# Write results to file (if applicable)
if (args$output_list) {
    write.table(enrich,
                gsub(".png", ".txt", args$output),
                sep       = "\t",
                row.names = FALSE)
}
message("Done.")
