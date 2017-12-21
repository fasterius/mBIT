#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "devtools", "ggbiplot", "caret")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  message("Installing missing packages ...")
  tryCatch (silent = TRUE,
            install.packages(setdiff(packages, rownames(installed.packages())),
                             repos = "http://cran.us.r-project.org"),
            warning = function(bc) {
              suppressPackageStartupMessages(library("devtools"))
              install_github("vqv/ggbiplot")
            },
            error = function(bc) {
              suppressPackageStartupMessages(library("devtools"))
              install_github("vqv/ggbiplot")
            })
}

# Command parser
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = paste("The columns of the input file should",
    "file should contain the names of the samples to be analysed on the form", 
    "<samplename_[abc]> in the columns. If samples are not specified, all of",
    "the samples in the file will be used."))
parser$add_argument("input",
                    type = "character",
                    help = "input file path")
parser$add_argument("output",
                    type = "character",
                    help = "output file path")
parser$add_argument("-l", "--log",
                    action = "store_true",
                    dest = "log",
                    help = "use log2(data) for analysis")
parser$add_argument("-H", "--hide-legend",
                    action = "store_false",
                    dest = "legend",
                    help = "hide legend")
parser$add_argument("-S", "--summary",
                    action = "store_true",
                    dest = "summary",
                    help = "show PCA summary")
parser$add_argument("-t", "--text",
                    action  = "store_true",
                    dest    = "text",
                    help    = "show text labels")
parser$add_argument("-c", "--coding",
                    type = "character",
                    dest = "coding",
                    default = "",
                    metavar = "",
                    help = "only include coding genes [biomaRt info file]")
parser$add_argument("-s", "--samples",
                    type = "character",
                    default = "all", 
                    dest = "input_samples",
                    metavar = "",
                    help = "sample names [format: <s1,s2,...,sN>]")
parser$add_argument("-p", "--palette",
                    type = "character",
                    dest = "palette",
                    default = paste0("#4e8ce4,#4ee4e4,#4ee499,#e4bf4e,",
                                     "#e4734e,#e44e4e,#e44e99,#734ee4,",
                                     "#999999"),
                    metavar = "",
                    help = "colour palette [format: <col1,col2,...,colN>]")
parser$add_argument("-f", "--filter",
                    type = "integer",
                    dest = "filter", 
                    default = 0,
                    metavar = "",
                    help = "filter data below threshold (integer)")
parser$add_argument("-C", "--components",
                    type = "character",
                    default = "1,2",
                    dest = "components",
                    metavar = "",
                    help = "the principal compontents to plot [default: 1,2]")
parser$add_argument("-g", "--groups-from-file",
                    type = "character",
                    dest = "groups_from_file", 
                    default = "",
                    metavar = "",
                    help = "2-column file containing sample groups")
parser$add_argument("-i", "--input-2",
                    type = "character",
                    dest = "input_2",
                    default = "none",
                    metavar = "",
                    help = "second input file")
args <- parser$parse_args()

# Functions -------------------------------------------------------------------

# Function for optional filtering of data
filter_data <- function(data) {
    
    # Filter for desired samples (if applicable)
    if (args$input_samples != "all") {
        samples <- gsub(",", "|", args$input_samples)
        data <- data[, grep(samples, names(data), value = TRUE)]
    }

    # Filter for coding genes (if applicable)
    if (args$coding != "") {
        # Read biomaRt info file
        message("Reading biomaRt data ...")
        info <- read.table(args$coding,
                           sep              = "\t",
                           header           = TRUE,
                           stringsAsFactors = FALSE)
        info <- info[c("ensembl_transcript_id",
                       "ensembl_gene_id",
                       "gene_biotype")]

        # Filter genes
        message("Removing non-coding genes ...")
        coding <- info[info$gene_biotype == "protein_coding", ]
        data <- data[(row.names(data) %in% coding$ensembl_transcript_id) |
                     (row.names(data) %in% coding$ensembl_gene_id), ]
    }

    # Filter on expression levels (if applicable)
    if (args$filter != 0) {
        start.genes <- length(data$ENSGID)
        data$mean <- rowMeans(subset(data, select = -ENSGID))
        data <- subset(data, mean > args$filter, select = -mean)
        end.genes <- length(data$ENSGID)
        cat(paste("Removed", start.genes - end.genes, 
                  "genes below cutoff\n"))
    }

    # Return filtered data
    return(data)
}

# Function for performing PCA
pc_analysis <- function(data) {

    # Transpose and convert to data.frame
    data_t <- as.data.frame(t(data))

    # Remove colums with no variance
    message("Removing zero-variance data ...")
    suppressPackageStartupMessages(suppressWarnings(library("caret")))
    zero_var <- nearZeroVar(data_t, names = TRUE)
    data_t2 <- data_t[, !(names(data_t) %in% zero_var)]
    message(paste("Removed", dim(data_t)[2] - dim(data_t2)[2], "/",
                  dim(data_t)[2], "columns without variance from the model."))

    # PCA analysis
    message("Computing principal components ...")
    pca <- prcomp(data_t2, center = TRUE, scale. = TRUE)

    # Return PCA
    return(pca)
}

# Function for plotting PCAs
plot_pca <- function(data, pca) {

    # Choice of components
    comps <- strtoi(strsplit(args$components, ",")[[1]])

    # Proportion of variance explained (for axis labels)
    pve <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
    pve_1 <- pve[comps[1]]
    pve_2 <- pve[comps[2]]

    # Set labels
    cols <- names(data)
    labels <- substr(cols, nchar(cols), nchar(cols))

    # Set groups
    if (args$groups_from_file != "") {

        # Groups from file
        groups <- read.table(args$groups, header = TRUE, sep = "\t")
        groups <- groups[, 2]
    } else {

        # Groups from sample names
        groups <- substr(cols, 1, nchar(cols) - 2)
    }

    # Set group order
    if (args$input_samples != "all") {
        inputs <- strsplit(args$input_samples, ",")[[1]]
    } else {
        inputs <- unique(groups)
    }
    
    # Factorise groups
    groups <- factor(groups, levels = inputs)
    
    # Axis labels
    lab_x <- paste0("PC", comps[1], " (", pve_1, " % of variance explained)")
    lab_y <- paste0("PC", comps[2], " (", pve_2, " % of variance explained)")

    # PCA plot
    suppressPackageStartupMessages(suppressWarnings(library("ggbiplot")))
    gg <- ggbiplot(pca,
                   obs.scale = 1,
                   var.axes  = FALSE,
                   groups    = groups,
                   choices   = c(comps[1], comps[2])) +
        geom_point(aes(colour  = groups), size = 3.5) +
        labs(x = lab_x, y = lab_y) +
        theme_bw()
    
    # Text
    if (args$text) {
        gg <- gg + geom_text(label  = labels,
                             colour = "white",
                             size   = 2.75)
    }

    # Colours
    palette <- strsplit(args$palette, ",")[[1]]

    if (length(unique(groups)) <= length(palette)) {
        gg <- gg +
            scale_colour_manual(name   = "",
                                values = palette)
    } else {
        message("Too few colours; using the default ggplot palette.")
        gg <- gg + labs(colour = "")
    }

    # Legend
    if (!args$legend) {
        gg <- gg + theme(legend.position = "none")
    }

    # Return graphical object
    return(gg)
}

# Analysis --------------------------------------------------------------------

# Load data
message("Reading data ...")
data <- read.table(args$input,
                   sep              = "\t",
                   row.names        = 1,
                   header           = TRUE,
                   stringsAsFactors = FALSE)

# Load and merge with secondary input data (if applicable)
if (args$input_2 != "none") {
    data_2 <- read.table(args$input_2,
                         header           = TRUE,
                         sep              = "\t",
                         row.names        = 1,
                         stringsAsFactors = FALSE)
    data <- merge(data, data_2, by = "row.names")
    row.names(data) <- data$Row.names
    data$Row.names <- NULL
}

# Filter data (if applicable)
data <- filter_data(data)

# Convert to numeric
for (col in names(data)) {
    data[col] <- suppressWarnings(as.numeric(data[[col]]))
}

# Use log(values) (if applicable)
if (args$log) {
    data <- log2(data + 1)
}

# Principal component analysis
pca <- pc_analysis(data)

# Show PCA summary (if applicable)
if (args$summary) {
    summary(pca)
}

# Plot pca
gg <- plot_pca(data, pca)

# Save to file
ggsave(args$output, gg, dpi = 300, width = 7, height = 5)
message("Done.")
