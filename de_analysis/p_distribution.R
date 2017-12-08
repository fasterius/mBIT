#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "ggplot2")
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

# Command parser
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = "Plot p-value distribution of DEG file")
parser$add_argument("input", type = "character", help = "input file path")
parser$add_argument("output", type = "character", help = "output file path")
args <- parser$parse_args()

# Read data
message("Reading data ...")
data <- read.table(args$input,
                   header           = TRUE,
                   sep              = "\t",
                   row.names        = 1, 
                   stringsAsFactors = FALSE)

# Plot histogram
message("Plotting data ...")
suppressPackageStartupMessages(suppressWarnings(library("ggplot2")))
gg <- ggplot(data = data, aes(p_value)) +
      geom_histogram(bins = 50, fill = "#4e8ce4", col = "black") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(x     = "P-values",
           y     = "Count",
           title = paste0("P-value distribution (", 
                         gsub(".txt", "", args$input), ")")) 

# Save and open image
ggsave(args$output, gg, dpi = 300, width = 7, height = 7)
message("Done.")
