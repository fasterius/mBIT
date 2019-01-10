#!/usr/local/bin/Rscript

# Install missing packages (if applicable)
packages <- c("argparse")
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

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = "Calculates sample RNA expression means.")
parser$add_argument("input",
                    type = "character",
                    help = "input data file path")
parser$add_argument("samples",
                    type = "character", 
                    help = "samples [format: <s1,s2,...,sN>]")
parser$add_argument("output",
                    type = "character",
                    help = "output data file path")
args <- parser$parse_args()

# Read input data
message("Reading data ...")
data <- read.table(args$input,
                   sep              = "\t",
                   header           = TRUE,
                   row.names        = 1,
                   stringsAsFactors = FALSE)

# Check for feature type
if (grepl("ENST", row.names(data)[1])) {
    id_name <- "ENSTID"
} else {
    id_name <- "ENSGID"
}

# Convert to numeric
message("Calculating sample means ...")
for (col in names(data)) {
    data[col] <- as.numeric(data[[col]])
}

# Calculate sample means
samples <- strsplit(args$samples, ",")[[1]]
for (sample in samples) {
    data[sample] <- rowMeans(data[c(grep(sample, names(data)))])
}

# Remove non-mean columns
data <- data[samples]
data[id_name] <- row.names(data)
data <- data[c(id_name, samples)]

# Remove NA-only rows (i.e. last row)
data <- data[rowSums(is.na(data)) != ncol(data), ]

# Save to file
write.table(data,
            file      = args$output,
            sep       = "\t",
            row.names = FALSE)
message("Done.")
