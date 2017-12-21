#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    message("Installing missing packages ...")
    install.packages(setdiff(packages, rownames(installed.packages())),
                 repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = "Add metadata to long-format data.")
parser$add_argument("input",
                    type    = "character",
                    help    = "input data file path")
parser$add_argument("metadata",
                    type    = "character",
                    help    = "metadata file path")
parser$add_argument("output",
                    type    = "character",
                    help    = "output data file path")
parser$add_argument("-m", "--merge-column",
                    type    = "character",
                    dest    = "merge_column",
                    default = "sample",
                    metavar = "",
                    help    = "metadata column to merge on [default: sample]")
parser$add_argument("-M", "--meta-columns",
                    type    = "character",
                    dest    = "meta_columns",
                    default = "all",
                    metavar = "",
                    help    = "columns to include [format: <col1,col2,...>")
args <- parser$parse_args()

# Read input data
message("Reading data ...")
data <- read.table(args$input,
                   sep              = "\t",
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Read metdata
message("Reading metadata ...")
metadata <- read.table(args$metadata,
                       sep              = "\t",
                       header           = TRUE,
                       stringsAsFactors = FALSE)

# Remove unwanted metadata columns (if applicable)
if (args$meta_columns != "all") {

    # Columns to keep
    to_keep <- strsplit(args$meta_columns, ",")[[1]]

    # Remove unwated columns
    metadata <- metadata[c(args$merge_column, to_keep)]
}

# Merge sample 1
message("Merging data and metadata ...")
merged <- merge(data, metadata,
                by.x = "sample_1",
                by.y = args$merge_column)

# Merge sample 2
merged <- merge(merged, metadata,
                by.x = "sample_2",
                by.y = args$merge_column)

# Fix column names
meta_cols <- gsub("\\.x", "", grep("\\.x", names(merged), value = TRUE))
meta_cols <- paste(meta_cols,
                   rep(c("2", "1"), each = length(meta_cols)),
                   sep = "_")
names(merged) <- c(names(data), meta_cols)

# Column ordering
merged <- merged[c("sample_1", "sample_2", sort(meta_cols), "r2")]

# Save to file
write.table(merged,
            args$output,
            sep = "\t",
            row.names = FALSE)
message("Done.")
