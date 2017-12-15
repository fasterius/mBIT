#!/usr/bin/env Rscript

# Command parser and packages -------------------------------------------------

# Install missing packages (if applicable)
packages <- c("argparse", "pheatmap")
if ( length(setdiff(packages, rownames(installed.packages()))) > 0 ) {
  cat("installing missing packages ...\n")
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
parser <- ArgumentParser(epilog = paste("Creates a correlation heatmap from",
                                        "expression estimates of RNA-seq", 
                                        "data."))
parser$add_argument("input",
                    type    = "character",
                    help    = "input data file path")
parser$add_argument("output",
                    type    = "character",
                    help    = "output data file path")
parser$add_argument("-a", "--annotate",
                    action  = "store_true",
                    dest    = "annotate",
                    help    = "annotate tiles with correlation coefficients")
parser$add_argument("-c", "--cluster",
                    action  = "store_true",
                    dest    = "cluster",
                    help    = "cluster the data")
parser$add_argument("-m", "--method",
                    type    = "character",
                    dest    = "method",
                    default = "pearson",
                    metavar = "",
                    help    = "correlation method [default: pearson]")
args <- parser$parse_args()

# Read expression data
message("Reading data ...")
data <- read.table(args$input,
                   row.names        = 1,
                   sep              = "\t",
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Convert to numeric
for (col in names(data)) {
    data[col] <- as.numeric(data[[col]])
}

# Logarithmise
data <- log2(data + 1)

# Calculate correlations
message("Calculating correlations ...")
correlations <- data.frame(row.names = names(data))
for (col_1 in names(data)) {

    col_cor <- vector()

    for (col_2 in names(data)) {
        vector_cor <- cor(x      = data[[col_1]],
                          y      = data[[col_2]],
                          method = args$method, 
                          use    = "complete.obs")
        col_cor <- c(col_cor, vector_cor)
    }

    correlations[col_1] <- col_cor

}

# Convert to matrix
correlations <- as.matrix(correlations)

# Plot the heatmap
message("Plotting ...")
suppressPackageStartupMessages(library("pheatmap"))
pheatmap(correlations,
         cluster_rows             = args$cluster,
         cluster_cols             = args$cluster,
         display_number           = args$annotate,
         clustering_distance_rows = "correlation",
         filename                 = args$output)

# Final message
message("Done.")
