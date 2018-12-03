#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    message("Installing missing packages ...")
    install.packages(setdiff(packages, rownames(installed.packages())),
                 repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = "Create long-format pairwise correlations.")
parser$add_argument("input",
                    type    = "character",
                    help    = "input data file path")
parser$add_argument("output",
                    type    = "character",
                    help    = "output data file path")
parser$add_argument("-l", "--log-normalise",
                    action  = "store_true",
                    dest    = "log",
                    help    = "log2-normalise data")
parser$add_argument("-m", "--method",
                    type    = "character",
                    dest    = "method",
                    default = "pearson",
                    metavar = "",
                    help    = "Pearson [default] or Spearman correlations")
args <- parser$parse_args()

# Functions -------------------------------------------------------------------

# Function for performing pairwise correlations
perform_correlations <- function(data) {

    # Initialise results dataframe
    corrs <- data.frame(sample_1         = character(),
                        sample_2         = character(),
                        r2               = numeric(),
                        stringsAsFactors = FALSE)

    # All samples to correlate
    samples <- names(data)

    # Total number of correlations to perform
    n_tot <- length(samples) * (length(samples) - 1) / 2 + length(samples)

    # Perform correlations
    for (s1 in samples) {
        for (s2 in samples) {

            # Check if mirrored comparison already exists
            s1s2 <- corrs[corrs$sample_1 == s1 & corrs$sample_2 == s2, ]
            s2s1 <- corrs[corrs$sample_1 == s2 & corrs$sample_2 == s1, ]
            if (nrow(s1s2) > 0 | nrow(s2s1) > 0) {
                next
            }

            # Progress message
            n <- nrow(corrs) + 1
            message("[", n, "/", n_tot, "]: Correlating ", s1, " and ", s2)

            # Correlate current samples
            if (s1 == s2) {
                r2 <- 1
            } else {
                corr <- cor(x      = as.numeric(data[[s1]]),
                            y      = as.numeric(data[[s2]]),
                            method = args$method)
                r2 <- corr * corr
            }

            # Add to correlations dataframe
            corrs[n, "sample_1"] <- s1
            corrs[n, "sample_2"] <- s2
            corrs[n, "r2"] <- r2
        }
    }

    # Return correlations
    return(corrs)
}

# Analysis --------------------------------------------------------------------

# Read input data
message("Reading data ...")
data <- read.table(args$input,
                   header           = TRUE,
                   sep              = "\t",
                   stringsAsFactors = FALSE,
                   row.names        = 1)

# Remove NA columns
data[sapply(data, function(x) all(is.na(x)))] <- NULL

# Log-normalise (if applicable)
if (args$log) {
    data <- log2(data + 1)
}

# Perform correlations
results <- perform_correlations(data)

# Save to file
write.table(results,
            args$output,
            sep = "\t",
            row.names = FALSE)
message("Done.")
