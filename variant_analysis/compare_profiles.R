#!/usr/bin/env Rscript

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
suppressPackageStartupMessages(library("seqCAT"))

# Get all SNV profiles in input
files <- list.files(args$input, full.names = TRUE)
profiles <- grep("profile.txt", files, value = TRUE)

# Read all profiles
profile_list <- list()
for (profile in profiles) {

    # Get current sample
    current_sample <- basename(profile)
    len <- length(current_sample)
    current_sample <- strsplit(current_sample[len], "\\.")[[1]][1]

    # Read current profile
    current_profile <- read_profile(profile, current_sample)

    # Append profile to profile list
    profile_list[[length(profile_list) + 1]] <- current_profile
}

# Compare all read profiles to each other
comparisons <- compare_many(profile_list)
similarities <- comparisons[[1]]

# Save similarities
write.table(similarities, args$output, sep = "\t", row.names = FALSE)
