#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "ggplot2")
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
parser <- ArgumentParser(epilog = "Create heatmaps from long-format data")
parser$add_argument("input",
                    type    = "character",
                    help    = "input data file path")
parser$add_argument("column",
                    type    = "character",
                    help    = "similarity measure column")
parser$add_argument("output",
                    type    = "character",
                    help    = "output image file path")
parser$add_argument("-H", "--hide-text",
                    action  = "store_true",
                    dest    = "hide_text",
                    help    = "hide text from plot")
parser$add_argument("-l", "--hide-legend",
                    action  = "store_true", 
                    dest    = "hide_legend",
                    help    = "hide legend from plot")
parser$add_argument("-c", "--cluster",
                    type    = "character",
                    dest    = "cluster",
                    default = "",
                    metavar = "",
                    help    = "cluster the resulting heatmap [all / group]")
parser$add_argument("-C", "--colour",
                    type    = "character",
                    dest    = "colour",
                    default = "#1954A6",
                    metavar = "",
                    help    = "colour base for gradient (default: #1954A6)")
parser$add_argument("-g", "--groups",
                    type    = "character",
                    dest    = "groups",
                    default = "",
                    metavar = "",
                    help    = "column containing sample groups")
parser$add_argument("-L", "--limits",
                    type    = "character",
                    dest    = "limits", 
                    metavar = "",
                    default = "0,0.5,0.9,1",
                    help    = "gradient limits")
parser$add_argument("-f", "--filter",
                    type    = "integer",
                    dest    = "filter",
                    default = 1000,
                    metavar = "",
                    help    = "max number of 0-columns")
parser$add_argument("-p", "--prefix",
                    type    = "character",
                    dest    = "prefix",
                    default = "sample",
                    metavar = "",
                    help    = "sample name prefix [default: sample]")
parser$add_argument("-s", "--size",
                    type    = "character",
                    dest    = "size",
                    default = "7x7",
                    metavar = "",
                    help    = "image height x width [default: 7x7]")
parser$add_argument("-t", "--title",
                    type    = "character",
                    dest    = "title",
                    default = "",
                    metavar = "",
                    help    = "plot title")
parser$add_argument("-z", "--legend-size",
                    type    = "double",
                    dest    = "legend_size",
                    default = 1,
                    metavar = "",
                    help    = "legend size (in inches)")
args <- parser$parse_args()

# Functions -------------------------------------------------------------------

# Function for filtering data
filter_data <- function(data) {

    # Remove samples below threshold (if applicable)
    if (args$filter != "") {

        # For every sample ...
        for (sample in unique(data$sample_1)) {

            # Calculate number of zero overlapping calls
            current <- data[data$sample_1 == sample, ]
            zeroes <- current[current$calls <= 10, ]
            
            # Remove sample from data if number of zeroes is above threshold
            if (nrow(zeroes) > args$filter) {
                data <- data[data$sample_1 != sample & 
                             data$sample_2 != sample, ]
            }
        }
    }

    # Return filtered data
    return(data)
}

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
    return(df_mat)
}

# Function for clustering data
cluster_data <- function(data) {

    # Cluster [all] or per [group]
    if (args$cluster == "all") {

        # Cluster data
        df_mat <- create_distance_matrix(data)
        clust <- hclust(dist(df_mat))
        order <- names(df_mat[clust$order])

    } else if (args$cluster == "group") {

        # Initialise order list
        order <- list()

        # Group column
        group_col <- paste0(args$groups, "_1")

        # Loop through each group
        for (group in unique(data[[group_col]])) {

            # Subset for current group
            current <- data[data[[group_col]] == group, ]

            # Cluster
            df_mat <- create_distance_matrix(current)
            clust <- hclust(dist(df_mat))
            order_current <- names(df_mat[clust$order])
            order <- c(order, order_current)
        }
    }

    # New cluster-based order
    data$s1 <- factor(data$s1, levels = order)
    data$s2 <- factor(data$s2, levels = order)

    # Return clustered data
    return(data)
}

# Function to mirror data to achieve a full square heatmap
mirror_data <- function(data) {

    # Loop through each sample combination
    for (s1 in unique(data$s1)) {
        for (s2 in unique(data$s1)) {

            # Add mirrored data where missing
            current <- data[data$s1 == s2 & data$s2 == s1, ]
            if (nrow(current) == 0) {
                data[nrow(data) + 1, "s1"] <- s2
                data[nrow(data), "s2"] <- s1
                data[nrow(data), args$column] <- data[data$s1 == s1 &
                                                      data$s2 == s2,
                                                          args$column]
            }
        }
    }

    # Return data
    return(data)
}

# Function for plotting the heatmap
plot_heatmap <- function(data) {

    # Get colour gradient limits
    lims <- as.numeric(strsplit(args$limits, ",")[[1]])
    limits <- c(lims[1], lims[2], lims[3], lims[4])

    # Gradient name
    gradient_name <- gsub("_", " ", args$column)
    gradient_name <- paste0(toupper(substr(gradient_name, 1, 1)),
                            substr(gradient_name, 2, nchar(gradient_name)))

    # Plotting
    gg <- ggplot(data, aes_string(x    = "s1",
                                  y    = "s2",
                                  fill = args$column)) + 
        geom_tile(colour = "white", size = 0.3) + 
        coord_equal() + 
        theme(axis.ticks       = element_blank(), 
              panel.background = element_blank(),
              axis.text.x      = element_text(angle = 90,
                                              hjust = 1,
                                              vjust = 0.5)) +
        labs(x = NULL,
             y = NULL,
             fill = gradient_name) + 
        scale_fill_gradientn(
            colours = c("white", "white", "#808080", args$colour), 
            limits  = c(0, max(limits)),
            values  = rescale(limits))

    # Add text (if applicable)
    if (!args$hide_text) {
        data[args$column] <- round(data[[args$column]], 2)
        gg <- gg + geom_text(data   = data,
                             colour = "white",
                             size   = 2,
                             aes_string(label = args$column))
    }

    # Hide legend (if applicable)
    if (args$hide_legend) {
        gg <- gg + theme(legend.position = "none")
    }

    # Change legend_size (if applicable)
    if (args$legend_size != "") {
        gg <- gg + theme(legend.key.height = unit(args$legend_size, "in"))
    }
      
    # Plot title (if applicable)
    if (args$title != "") {
        gg <- gg + ggtitle(args$title) +
            theme(plot.title = element_text(hjust = 0.5))
    }

    # Return graphical object
    return(gg)
}

# Analysis --------------------------------------------------------------------

# Read data
message("Reading data ...")
data <- read.table(args$input,
                   sep              = "\t",
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Load packages
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("dplyr"))

# Get sample columns
data$s1 <- data[, paste0(args$prefix, "_1")]
data$s2 <- data[, paste0(args$prefix, "_2")]

# Filter data
data <- filter_data(data)

# Cluster data (if applicable)
if (args$cluster != "") {
    message("Clustering data ...")
    data <- cluster_data(data)
}

# More data to get full square heatmap
message("Plotting heatmap ...")
data <- mirror_data(data)

# Plot
gg <- plot_heatmap(data)

# Save to file
size <- as.numeric(strsplit(args$size, "x")[[1]])
ggsave(args$output, gg, dpi = 300, height = size[1], width = size[2])
message("Done.")
