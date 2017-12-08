#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "ggplot2")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  message("Installing missing packages ...")
  tryCatch (silent = true,
            install.packages(setdiff(packages, rownames(installed.packages())),
                             repos = "http://cran.us.r-project.org"),
            warning = function(bc) {
              source("http://bioconductor.org/bioclite.r")
              bioclite(setdiff(packages, rownames(installed.packages())))
            },
            error = function(bc) {
              source("http://bioconductor.org/bioclite.r")
              bioclite(setdiff(packages, rownames(installed.packages())))
            })
}

# Command parser
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = gsub("[\r\n]", "", ""))
parser$add_argument("input",
                    type    = "character",
                    help    = "input DEG file path")
parser$add_argument("output",
                    type    = "character",
                    help    = "output image file path")
parser$add_argument("-f", "--FDR-cutoff",
                    type    = "double",
                    dest    = "fdr_cutoff",
                    default = 0.01,
                    metavar = "",
                    help    = "FDR cutoff for significance [default: 0.01]")
parser$add_argument("-F", "--FC-cutoff",
                    type    = "integer",
                    dest    = "fc_cutoff",
                    default = 2,
                    metavar = "",
                    help    = "FC cutoff for significance [default: 2]")
parser$add_argument("-c", "--colour",
                    type    = "character",
                    dest    = "colour", 
                    default = "#e4bf4e,#4e8ce4",
                    metavar = "", 
                    help    = "colours [format: <up,down>]")
args <- parser$parse_args()

# Read data
message("Reading data ...")
data <- read.table(args$input,
                   row.names        = 1,
                   sep              = "\t",
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Add -log10(FDR) to data
data$logFDR <- -log10(data$FDR)
data <- na.omit(data[c("log2FC", "logFDR")])

# Set significance- and FC-dependent colours for each point
data$group <- "non_significant"
# data[data$logFDR < -log10(args$fdr_cutoff), "group"] <- "non_significant"
data[data$logFDR >= -log10(args$fdr_cutoff) &
     data$log2FC >= log2(args$fc_cutoff), "group"] <- "up"
data[data$logFDR >= -log10(args$fdr_cutoff) &
     data$log2FC <= -log2(args$fc_cutoff), "group"] <- "down"
data$group <- factor(data$group, levels = c("up", "down", "non_significant"))

# Colours
message("Plotting data ...")
suppressPackageStartupMessages(library("ggplot2"))
colours <- strsplit(args$colour, ",")[[1]]

# Plot
gg <- ggplot(data, aes(x = log2FC, y = logFDR, colour = group)) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = log2(args$fc_cutoff),
               linetype   = "dotted",
               size       = 0.5,
               colour     = "#999999") +
    geom_vline(xintercept = -log2(args$fc_cutoff),
               linetype   = "dotted",
               size       = 0.5,
               colour     = "#999999") +
    geom_hline(yintercept = -log10(args$fdr_cutoff),
               linetype   = "dotted",
               size       = 0.5,
               colour     = "#999999") +
    theme_bw() +
    labs(title = paste0("Volcano plot (", args$input, ")"),
         x     = "log2(Fold Change)",
         y     = "-log10(FDR)") +
    theme(plot.title      = element_text(hjust = 0.5),
          legend.position = "none") +
    scale_colour_manual(values = c(colours[1], colours[2], "#cccccc"))

# Save plot
ggsave(args$output, gg, dpi = 300, width = 7, height = 7)
message("Done.")
