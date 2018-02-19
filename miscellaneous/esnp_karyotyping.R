#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "tidyr", "dplyr", "eSNPKaryotyping", "seqCAT",
              "devtools")
if ( length(setdiff(packages, rownames(installed.packages()))) > 0 ) {
    message("Installing missing packages ...")
    tryCatch (silent = TRUE,
            install.packages(setdiff(packages, rownames(installed.packages())),
                             repos = "http://cran.us.r-project.org"),
            warning = function(bc) {
                source("http://bioconductor.org/biocLite.R")
                biocLite(setdiff(packages, rownames(installed.packages())))
                devtools::install_github(
                    "BenvenLab/eSNPKaryotyping/eSNPKaryotyping")
            },
            error = function(bc) {
                source("http://bioconductor.org/biocLite.R")
                biocLite(setdiff(packages, rownames(installed.packages())))
                devtools::install_github(
                    "BenvenLab/eSNPKaryotyping/eSNPKaryotyping")
            })
}

# Command parser
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = paste("Perform eSNP-Karyotyping on SNV",
                                        "profiles."))
parser$add_argument("input",
                    type    = "character",
                    help    = "input data file or directory path")
parser$add_argument("output",
                    type    = "character",
                    help    = "output file path for results table")
parser$add_argument("-d", "--directory",
                    action  = "store_true",
                    dest    = "directory",
                    help    = "input is a directory with several SNV profiles")
args <- parser$parse_args()

# Functions -------------------------------------------------------------------

# Function for reading one or more SNV profiles
read_profiles <- function(input, dir) {

    # Initiate profile list
    profiles <- list()

    # List profiles
    if (!dir) {
        
        # Single input
        files <- list.files(dirname(input), pattern = input)

    } else {
        
        # Directory input
        files <- list.files(input, pattern = "\\.profile\\.txt")
    }

    # Read profile(s)
    nn_total = length(files)
    nn = 1
    for (file in files) {

        # Get sample
        sample <- gsub("\\.profile\\.txt", "", basename(file))

        # Read profile
        message(paste0("Reading profile for ", sample, " in file ", file, " [",
                       nn, " / ", nn_total, "]"))
        profile <- suppressMessages(as.data.frame(read_profile(file, sample)))

        # Cleanup
        profile <- profile[
            c("sample", "seqnames", "start", "DP", "AD1", "AD2")]
        names(profile) <- c("sample", "chr", "position", "DP", "AD1", "AD2")
        profile <- unique(profile)

        # Add profile to list
        profiles[[length(profiles) + 1]] <- profile

        # Increment counter
        nn = nn + 1
    }

    # Return all profiles
    return(profiles)
}

# Function for performing eSNP-Karyotyping
esnp_karyotyping <- function(profiles) {

    # Initiate results dataframe
    results <- data.frame(sample           = character(),
                          total_SNVs       = character(),
                          mean_AR          = numeric(),
                          chr1             = numeric(),
                          chr2             = numeric(),
                          chr3             = numeric(),
                          chr4             = numeric(),
                          chr5             = numeric(),
                          chr6             = numeric(),
                          chr7             = numeric(),
                          chr8             = numeric(),
                          chr9             = numeric(),
                          chr10            = numeric(),
                          chr11            = numeric(),
                          chr12            = numeric(),
                          chr13            = numeric(),
                          chr14            = numeric(),
                          chr15            = numeric(),
                          chr16            = numeric(),
                          chr17            = numeric(),
                          chr18            = numeric(),
                          chr19            = numeric(),
                          chr20            = numeric(),
                          chr21            = numeric(),
                          chr22            = numeric(),
                          chrX             = numeric(),
                          chrY             = numeric(),
                          stringsAsFactors = FALSE)

    # Loop through each profile
    for (profile in profiles) {
        
        # Get sample
        sample <- unique(profile$sample)
        profile$sample <- NULL

        # Perform eSNP-Karyotyping (globally)
        major_minor <- MajorMinorCalc(profile, 10, 1000000, 0.2)
        total_snvs <- nrow(major_minor)
        mean_ar <- mean(major_minor$MajorMinor)
        
        # Per chromosome
        chr <- major_minor %>%
            group_by(chr) %>%
            summarise(mean_ar = mean(MajorMinor))
        chr <- as.data.frame(chr, stringsAsFactors = FALSE)

        # Add to results
        n <- nrow(results) + 1
        results[n, "sample"] <- sample
        results[n, "total_SNVs"] <- total_snvs 
        results[n, "mean_AR"] <- mean_ar 
        for (col in names(results)[4:27]) {

            # Get current chromosome data
            current <- chr[chr$chr == gsub("chr", "", col), "mean_ar"]
            
            # Check if data for current chromosome exists
            if (length(current) == 0) {
                current <- NA
            }

            # Add to results
            results[n, col] <- current
        }
    }

    # Return results
    return(results)
}

# Analysis --------------------------------------------------------------------

# Load packages
message("Loading packages ...")
suppressPackageStartupMessages(library("seqCAT"))
suppressPackageStartupMessages(library("eSNPKaryotyping"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))

# Read data
profiles <- read_profiles(args$input, args$directory)

# Perform eSNP-Karyotyping
message("Performing eSNP-Karyotyping ...")
results <- esnp_karyotyping(profiles)

# Write to file
write.table(results, args$output, sep = "\t", row.names = FALSE)
message("Done.")
