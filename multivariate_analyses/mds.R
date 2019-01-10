#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "ggplot2", "smacof")
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
parser <- ArgumentParser(epilog = "Create MDS from wide-format data.")
parser$add_argument("input",
                    type    = "character",
                    help    = "input data file path")
parser$add_argument("column",
                    type    = "character",
                    help    = "similarity measure column")
parser$add_argument("output",
                    type    = "character",
                    help    = "output image file path")
parser$add_argument("-l", "--hide-legend",
                    action  = "store_false",
                    dest    = "legend",
                    help    = "hide plot legend")
parser$add_argument("-f", "--filter",
                    type    = "integer",
                    dest    = "filter",
                    metavar = "",
                    default = 1000,
                    help    = "max number of 0-columns")
parser$add_argument("-g", "--groups",
                    type    = "character",
                    dest    = "groups",
                    default = "",
                    metavar = "",
                    help    = "known sample groups column")
parser$add_argument("-p", "--prefix",
                    type    = "character",
                    dest    = "prefix",
                    default = "sample",
                    metavar = "",
                    help    = "sample column prefix [default: sample]")
parser$add_argument("-P", "--palette",
                    type = "character",
                    dest = "palette",
                    default = paste0("#4e8ce4,#4ee4e4,#4ee499,#e4bf4e,",
                                     "#e4734e,#e44e4e,#e44e99,#734ee4,",
                                     "#999999"),
                    metavar = "",
                    help = "colour palette [format: <col1,col2,...,colN>]")
parser$add_argument("-s", "--size",
                    type    = "character",
                    dest    = "size",
                    metavar = "",
                    default = "7x7",
                    help    = "canvas height x width [default: 7x7]")
parser$add_argument("-m", "--method",
                    type    = "character",
                    dest    = "method",
                    default = "smacof",
                    metavar = "",
                    help    = "MDS method: smacof [default] or cmdscale")
parser$add_argument("-r", "--remove-samples",
                    type    = "character",
                    dest    = "remove_samples",
                    default = "",
                    metavar = "",
                    help    = "samples to remove [format: <s1,s2,...>]")
args <- parser$parse_args()

# Functions -------------------------------------------------------------------

# Function for filtering input data
filter_data <- function(data) {

    # Remove samples (if applicable)
    if (args$remove_samples != "") {
        remove_samples <- strsplit(args$remove_samples, ",")[[1]]
        data <- data[!(data$s1 %in% remove_samples), ]
        data <- data[!(data$s2 %in% remove_samples), ]
    }

    # Filter samples below threshold (if applicable)
    if (args$filter !=  1000) {

        # For every sample ...
        for (sample in unique(data$s1)) {

            # Calculate number of zero overlapping calls
            current <- data[data$s1 == sample, ]
            zeroes <- current[current$calls <=  10, ]

            # Remove sample from data if number of zeroes is above threshold
            if (nrow(zeroes) > args$filter) {
                data <- data[data$s1 !=  sample &
                             data$s2 !=  sample, ]
            }
        }
    }

    # Return filtered data
    return(data)
}

# Function for creating a distance matrix
create_distance_matrix <- function(data) {

    # Initialise empty matrix-like dataframe
    samples <- unique(data$s1)
    samples_n <- length(samples)
    df_mat <- as.data.frame(matrix(0, ncol = samples_n, nrow = samples_n))
    row.names(df_mat) <- samples
    names(df_mat) <- samples

    # Add data to the matrix-like data frame in a sample-per-sample manner
    for (s1 in samples) {
        for (s2 in samples) {
            current <- data[data$s1 == s1 & data$s2 == s2, args$column]
            if (length(current) == 0) {
                df_mat[s1, s2] <- data[data$s1 == s2 &
                                       data$s2 == s1, args$column]
            } else {
                df_mat[s1, s2] <- current
            }
        }
    }
    
    # Return distance matrix
    df_mat <- dist(df_mat)
    return(df_mat)
}

# Function for MDS
plot_mds <- function(data, df_mat) {

    # Multidimensional scaling
    if (args$method == "cmdscale") {

        # Use CMD-scale
        mds <- as.data.frame(cmdscale(df_mat))

    } else if (args$method == "smacof") {

        # Use SMACOF
        suppressPackageStartupMessages(library("smacof"))
        mds <- as.data.frame(smacofSym(df_mat, ndim <- 2)$conf)
    }

    # Add groups (if applicable)
    if (args$groups != "") {
        
        # Get group column
        group_col <- paste0(args$groups, "_1")
        data[[group_col]] <- factor(data[[group_col]])

        # Merge groups with MDS
        mds <- merge(mds, 
                     data[c("s1", group_col)],
                     by.x  = "row.names",
                     by.y  = "s1",
                     all.x = TRUE)
        mds <- unique(mds)
        names(mds) <- c("sample", "x", "y", "group")

    } else {

        # Add non-existing groups
        mds$group <- factor(1)
        mds$sample <- row.names(mds)
        names(mds) <- c("x", "y", "group", "sample")
    }

    # Plot
    suppressPackageStartupMessages(library("ggplot2"))
    gg <- ggplot(mds, aes(x      = x,
                          y      = y,
                          colour = group)) +
        geom_point(size = 3) +
        labs(x      = "Dimension 1",
             y      = "Dimension 2",
             colour = "",
             title  = paste0("Multidimensional scaling [",
                             args$method, "]")) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text  = element_blank(),
              axis.ticks = element_blank())

    # Colour palette
    palette <- strsplit(args$palette, ",")[[1]]
    if (length(unique(mds$group)) <= length(palette)) {
        gg <- gg +
            scale_colour_manual(values = palette)
    }

    # Hide legend (if applicable)
    if (!args$legend) {
        gg <- gg + theme(legend.position = "none")
    }

    # Return graphical object
    return(gg)
}

# Analysis --------------------------------------------------------------------

# Read data
data <- read.table(args$input,
                   sep              = "\t",
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Get sample columns
data$s1 <- data[, paste0(args$prefix, "_1")]
data$s2 <- data[, paste0(args$prefix, "_2")]

# Filter data (if applicable)
data <- filter_data(data)

# Create distance matrix
message("Clustering data ...")
df_mat <- create_distance_matrix(data)

# Create and plot MDS
message("Plotting ...")
gg <- plot_mds(data, df_mat)

# Save plot
size <- as.numeric(strsplit(args$size, "x")[[1]])
ggsave(args$output, gg, dpi = 300, height = size[1], width = size[2])
message("Done.")
