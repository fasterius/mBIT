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
parser <- ArgumentParser(epilog = "Compare SNV profiles in directory.")
parser$add_argument("input",
                    type    = "character",
                    help    = "input directory containing SNV profiles")
parser$add_argument("output",
                    type    = "character",
                    help    = "path for output similarities file")
args <- parser$parse_args()

# Load seqCAT
message("Loading packages ...")
suppressPackageStartupMessages(library("seqCAT"))

# Get all SNV profiles in input
files <- list.files(args$input, full.names = TRUE)
profiles <- grep("profile.txt", files, value = TRUE)

# Calculate total number of profiles to read and initialise counter
nn_tot <- length(profiles)
nn <- 1

# Read all profiles
profile_list <- list()
for (profile in profiles) {

    # Get current sample
    current_sample <- basename(profile)
    len <- length(current_sample)
    current_sample <- strsplit(current_sample[len], "\\.")[[1]][1]

    # Read current profile
    message(paste0("Reading profile for ", current_sample, " in file ",
                   basename(profile), " [", nn, " / ", nn_tot, "]"))
    current_profile <- suppressMessages(read_profile(profile, current_sample))

    # Append profile to profile list
    profile_list[[length(profile_list) + 1]] <- current_profile

    # Increment counter
    nn <- nn + 1
}

# Compare all read profiles to each other
comparisons <- compare_many(profile_list)
similarities <- comparisons[[1]]

# Save similarities
write.table(similarities, args$output, sep = "\t", row.names = FALSE)
