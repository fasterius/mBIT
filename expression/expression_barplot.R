#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("ggplot2", "dplyr", "reshape")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    message("Installing missing packages ...\n")
    tryCatch (silent = TRUE,
            install.packages(setdiff(packages, rownames(installed.packages())),
                             repos = "http://cran.us.r-project.org"))
}

# Command parser
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = paste("Plots the expression levels of the",
    "supplied gene(s) or isoform(s), using Ensembl IDs."))
parser$add_argument("input",
                    type    = "character",
                    help    = "input file path")
parser$add_argument("features",
                    type    = "character",
                    help    = "input features on the form <ID1,ID2,...,IDN>")
parser$add_argument("output",
                    type    = "character",
                    help    = "output file path")
parser$add_argument("-s", "--samples",
                    type    = "character",
                    dest    = "samples",
                    default = "", 
                    metavar = "",
                    help    = "samples on the form <s1,s2,...>")
parser$add_argument("-i", "--input-2",
                    type    = "character",
                    dest    = "input_2",
                    default = "",
                    metavar = "",
                    help    = "second input file")
args <- parser$parse_args()

# Function for converting data to numeric and excluding samples
fix_data <- function(data, features) {

    for (col in names(data)) {
        data[col] <- round(as.numeric(data[[col]]), 1)
    }

    # Only use specified samples (if applicable)
    if ( args$samples != "" ) {
        samples <- gsub(",", "|", args$samples)
        data <- data[, grep(samples, names(data), value = TRUE)]
    }

    # Select features
    data <- data[grep(features, row.names(data), value = TRUE), ]

    # Return data
    return(data)
}

# Define gene/isoform feature
features <- gsub(",", "|", args$features)
if (grepl("ENST", features)) {
    feature_type <- "Isoform"
} else if (grepl("ENSG", features)) {
    feature_type <- "Gene"
} else {
    stop("wrong feature type; please use ENSGID or ENSTID.")
}

# Load packages
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("reshape2"))

# Load data and make numeric
message("Reading data ...")
data <- read.table(args$input,
                   row.names        = 1, 
                   sep              = "\t",
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Fix data
data <- fix_data(data, features)

# Second input file (if applicable)
if (args$input_2 != "") {

    # Read second input
    data_2 <- read.table(args$input_2,
                       row.names        = 1, 
                       sep              = "\t",
                       header           = TRUE,
                       stringsAsFactors = FALSE)

    # Fix secondary data
    data_2 <- fix_data(data_2, features)

    # Merge with primary data
    data <- merge(data, data_2, by = "row.names")
    row.names(data) <- data$Row.names
    data$Row.names <- NULL
}

# Melt data
message("Plotting ...")
data[[feature_type]] <- row.names(data)
data <- melt(data, id = feature_type)
names(data) <- c(feature_type, "Sample", "TPM")

# Colour palette
palette <- c("#0d2d59", "#1954a6", "#4e8ce4", "#a6c6f2", # Blues
             "#59460d", "#a48119", "#e4bf4e", "#f2dfa6", # Yellows
             "#333333", "#666666", "#999999", "#cccccc") # Grays

# Plot isoform expression levels
gg <- ggplot(data, aes_string(x    = feature_type,
                              y    = "TPM",
                              fill = "Sample")) +
    geom_bar(stat     = "identity",
             position = "dodge",
             colour   = "black",
             size     = 0.3) +
    geom_text(aes(label = format(TPM, ndigits = 1)),
              position = position_dodge(width = 0.9),
              size     = 2.8,
              vjust    = -0.5) +
    theme_bw() +
    theme(axis.text   = element_text(size = 8),
          axis.title  = element_text(size = 12)) +
    labs(x = feature_type,
         y = "RNA expression (TPM)") +
    scale_fill_manual(values = palette, name = "")

# Save to file
ggsave(args$output, gg, dpi = 300, width = 7, height = 7)
message("Done.")
