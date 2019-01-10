#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "ggplot2", "pls")
if ( length(setdiff(packages, rownames(installed.packages()))) > 0 ) {
    message("Installing missing packages ...")
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
parser <- ArgumentParser(epilog = "Perform PLS on multivariate data.")
parser$add_argument("input",
                    type    = "character",
                    help    = "input data file path")
parser$add_argument("output",
                    type    = "character",
                    help    = "output image file path")
parser$add_argument("-p", "--print-coefs",
                    action  = "store_true",
                    dest    = "print_coefs",
                    help    = "print the variable coefficients")
parser$add_argument("-v", "--validation",
                    type    = "character",
                    dest    = "validation",
                    metavar = "",
                    default = "CV",
                    help    = "Validation strategy [CV]")
parser$add_argument("-s", "--size",
                    type    = "character",
                    dest    = "size",
                    metavar = "",
                    default = "7x7",
                    help = "canvas height x width [default: 7x7]")
args <- parser$parse_args()

# Functions -------------------------------------------------------------------

# Function for finding the ideal number of dimensions
find_ncomps <- function(data) {

    # Initialise the model
    pls_model <- plsr(group ~ ., data = data, validation = args$validation)

    # Find the ideal number of dimensions
    cv <- RMSEP(pls_model)
    ncomps <- which.min(cv$val[estimate = "adjCV", , ]) - 1

    # Return the ideal dimensions
    return(ncomps)
}

# Function for performing PLS
run_pls <- function(data, ncomps) {

    # Perform PLS with the specified number of dimensions
    pls_model <- plsr(group ~ ., data = data, ncomp = ncomps)

    # Get coefficients
    coefs <- coef(pls_model)
    sum_coefs <- sum(sapply(coefs, abs))
    coefs <- coefs / sum_coefs

    # Sort coefficients and convert to data frame
    coefs <- sort(coefs[, 1, 1])
    coefs <- data.frame(coefs)
    names(coefs) <- "coefficient"
    coefs$variable <- factor(row.names(coefs), levels = row.names(coefs))

    # Return dataframe with coefficients
    return(coefs)

}

# Function for plotting PLS coefficients
plot_pls <- function(coefs) {

    # Load ggplo2 package
    suppressPackageStartupMessages(library("ggplot2"))

    # Plot
    gg <- ggplot(coefs, aes(x = variable,
                            y = coefficient)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    labs(x = NULL,
         y = "Regression coefficient") +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5))

    # Return ggplot object
    return(gg)

}

# Analysis --------------------------------------------------------------------

# Read data
message("Reading data ...")
data <- read.table(args$input,
                   sep              = "\t",
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Set first column as "group"
names(data)[1] <- "group"

# Find ideal number of dimensions
message("Finding the ideal number of dimensions ...")
suppressPackageStartupMessages(library("pls"))
ncomps <- find_ncomps(data)

# Perform PLS and get coefficients
message("Performing PLS ...")
coefs <- run_pls(data, ncomps)

# Create a barplot
message("Plotting ...")
gg <- plot_pls(coefs)

# Save plot
size <- as.numeric(strsplit(args$size, "x")[[1]])
ggsave(args$output, gg, dpi = 300, height = size[1], width = size[2])

# Print coefficients (if applicable)
if (args$print_coefs) {

    # Add numbered row names
    row.names(coefs) <- NULL

    # Print coefficients and variables
    message("PLS coefficients:")
    print(coefs)
}
