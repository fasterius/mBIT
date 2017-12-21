#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "ggplot2", "ggdendro", "fossil")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
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
parser <- ArgumentParser(epilog = paste("Create dendrograms from long-format",
                                        " data."))
parser$add_argument("input",
                    type    = "character",
                    help    = "input data file path")
parser$add_argument("column",
                    type    = "character",
                    help    = "similarity measure column")
parser$add_argument("output",
                    type    = "character",
                    help    = "output image file path")
parser$add_argument("-f", "--filter",
                    type    = "integer", 
                    dest    = "filter",
                    default = 1000,
                    metavar = "",
                    help    = "max number of 0-columns [default: 1000]")
parser$add_argument("-g", "--groups",
                    type    = "character",
                    dest    = "groups",
                    default = "",
                    metavar = "",
                    help    = "known sample groups column")
parser$add_argument("-k", "--k-clusters",
                    type    = "integer",
                    dest    = "k",
                    default = 1,
                    metavar = "",
                    help    = "number of clusters [default: 1]")
parser$add_argument("-p", "--prefix",
                    type    = "character",
                    dest    = "prefix",
                    default = "sample",
                    metavar = "",
                    help    = "sample column name prefix [default: sample]")
parser$add_argument("-s", "--size",
                    type    = "character",
                    dest    = "size",
                    default = "7x7",
                    metavar = "",
                    help    = "canvas height x width [default: 7x7]")
parser$add_argument("-t", "--text-size",
                    type    = "double",
                    dest    = "text_size",
                    default = 4,
                    metavar = "",
                    help    = "text size [default: 4]")
parser$add_argument("-m", "--method",
                    type    = "character",
                    dest    = "method",
                    default = "complete",
                    metavar = "",
                    help    = paste("clustering method: complete, single,",
                                    "ward, centroid or average",
                                    "[default: complete]"))
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
    return(df_mat)
}

# Function for plotting dendrograms
plot_dendrogram <- function(data, clust, df_mat) {

    # Create contigency table for Adjusted Rand Index (if applicable)
    if (args$groups != "") {

        # True labels for ARI
        group_col <- paste0(args$groups, "_1")
        data[[group_col]] <- factor(data[[group_col]])
        true_labels <- as.numeric(unique(
            data[c(paste0(args$prefix, "_1"), group_col)])[[group_col]])

        # Calculate ARI
        suppressPackageStartupMessages(library("fossil"))
        ari <- adj.rand.index(true_labels, clust)
        ari <- format(ari, digits = 3)
        title <- paste0("Hierarchical clustering [", args$method,
                        "]\n(ARI = ", ari, ")")
    } else {

        # No ARI
        title <- paste0("Hierarchical clustering [", args$method, "]")
    }

    # Plot
    gg <- ggplot() +
        geom_segment(data = segment(dendr),
                     aes(x = x,
                         y = y,
                         xend = xend,
                         yend = yend)) +
        geom_text(data = label(dendr),
                  size = args$text_size,
                  aes(x, y,
                      label = label,
                      hjust = 0,
                      color = cluster)) +
        # geom_rect(data = rect, aes(xmin = X1 - 0.3, xmax = X2 + 0.3,
                                 # ymin = 0, ymax = ymax),
                # color = "#1954A6", fill = NA) +
        # geom_hline(yintercept = 0.33, color = "#cccccc") +
        coord_flip() +
        scale_y_reverse(expand = c(0.2, 0)) +
        theme_dendro() +
        theme(legend.position = "none",
              plot.title      = element_text(hjust = 0.5)) +
        labs(title = title)

    # Return graphical object
    return(gg)
}

# Analysis --------------------------------------------------------------------
 
# Load packages
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggdendro"))

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

# Cluster using preferred method
hca <- hclust(dist(df_mat), method = args$method)
k <- args$k
clust <- cutree(hca, k = k)

# Create dendrogram
dendr <- dendro_data(hca, type = "rectangle")
clust_df <- data.frame(label = row.names(df_mat), cluster = factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]], clust_df, by = "label")
rect <- aggregate(x ~ cluster, label(dendr), range)
rect <- data.frame(rect$cluster, rect$x)
ymax <- mean(hca$height[length(hca$height) - (k - 2):(k - 1)])

# Plot dendrogram
message("Plotting dendrogram ...")
gg <- plot_dendrogram(data, clust)

# Save to file
size <- as.numeric(strsplit(args$size, "x")[[1]])
ggsave(args$output, gg, dpi = 300, width = size[1], height = size[2])
message("Done.")
