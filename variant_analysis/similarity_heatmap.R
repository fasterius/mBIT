#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "seqCAT", "ggplot2")
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
parser <- ArgumentParser(epilog = "Create a heatmap of genetic similarities.")
parser$add_argument("input",
                    type    = "character",
                    help    = "input file path")
parser$add_argument("output",
                    type    = "character",
                    help    = "output image path")
parser$add_argument("-s", "--plot-size",
                    type    = "character",
                    dest    = "size",
                    metavar = "",
                    default = "7x7",
                    help    = "plot size [default: 7x7]")
parser$add_argument("-S", "--axis-size",
                    type    = "double",
                    dest    = "axis_size",
                    metavar = "",
                    default = 10,
                    help    = "axis text size [default: 10]")
parser$add_argument("-a", "--annotate",
                    action  = "store_false",
                    dest    = "annotate",
                    help    = "annotate cells [default: TRUE]")
parser$add_argument("-c", "--cluster",
                    action  = "store_true",
                    dest    = "cluster",
                    help    = "cluster the heatmap [default: FALSE]")
args <- parser$parse_args()

# Load packages
message("Loading packages ...")
suppressPackageStartupMessages(library("seqCAT"))
suppressPackageStartupMessages(library("ggplot2"))

# Read data
message("Reading data ...")
data <- read.table(args$input,
                   header           = TRUE,
                   sep              = "\t",
                   stringsAsFactors = FALSE)

# Plot heatmap
message("Plotting heatmap ...")
gg_heatmap <- plot_heatmap(data,
                           annotate = args$annotate,
                           cluster = args$cluster)
gg_heatmap <- gg_heatmap +
    theme(axis.text = element_text(size = args$axis_size))

# Save heatmap
size <- as.numeric(strsplit(args$size, "x")[[1]])
ggsave(args$output, gg_heatmap, dpi = 300, width = size[1], height = size[2])
message("Done.")
