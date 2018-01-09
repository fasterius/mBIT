#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "seqCAT")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
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

# Argument parser
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = paste0("Create SNV profiles for VCF files ",
                                         "in a directory. File name format ",
                                         "should be '<sample>.vcf[.gz]'."))
parser$add_argument("input",
                    type    = "character",
                    help    = "path to directory with VCF file(s)")
parser$add_argument("output",
                    type    = "character",
                    help    = "output path to put all SNV profiles")
parser$add_argument("-p", "--python",
                    action  = "store_true",
                    dest    = "python",
                    help    = "use the Python implementation instead of R")
parser$add_argument("-r", "--recursive",
                    action  = "store_true",
                    dest    = "recursive",
                    help    = "list all VCF files recursively")
parser$add_argument("-f", "--filter-depth",
                    type    = "integer",
                    dest    = "filter_depth",
                    metavar = "",
                    default = 10,
                    help    = "filter depth [default: 10]")
parser$add_argument("-P", "--pattern",
                    type    = "character",
                    dest    = "pattern",
                    metavar = "",
                    default = NULL,
                    help    = "pattern to filter files on (e.g. a single VCF)")
args <- parser$parse_args()

# Load seqCAT
suppressPackageStartupMessages(library("seqCAT"))

# List files to create profiles for
files <- list.files(path       = args$input,
                    full.names = TRUE,
                    pattern    = args$pattern,
                    recursive  = args$recursive)

# Get only VCF files
files <- grep(".vcf", files, value = TRUE)

# Create SNV profiles for each VCF
for (vcf in files) {

    # Get current sample
    sample <- strsplit(vcf, "/")[[1]]
    len <- length(sample)
    sample <- strsplit(sample[len], "\\.")[[1]][1]

    # Set current output
    output <- paste0(args$output, "/", sample, ".profile.txt")

    # Create profile for current input file
    tryCatch({
        create_profile(vcf,
                       sample,
                       output,
                       args$filter_depth,
                       args$python)
    }, error = function(e) {
        message("ERROR when creating profile for ", sample,
                "; continuing to the next sample.")
    })
}
