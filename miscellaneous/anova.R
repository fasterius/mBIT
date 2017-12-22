#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "tidyr", "ggplot2", "multcompView")
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
parser <- ArgumentParser(epilog = paste0("Perform ANOVA and plot both a",
                                         "boxplot and its confidence",
                                         "intervals"))
parser$add_argument("input",
                    type    = "character",
                    help    = "input data file path")
parser$add_argument("output",
                    type    = "character",
                    help    = "output image file path")
parser$add_argument("treatment",
                    type    = "character",
                    help    = "treatment variable")
parser$add_argument("response",
                    type    = "character",
                    help    = "response variable")
parser$add_argument("-c", "--conf-level",
                    type    = "double",
                    default = 0.95,
                    dest    = "conf_level",
                    metavar = "",
                    help    = "confidence level [default: 0.95]")
parser$add_argument("-p", "--palette",
                    type    = "character",
                    default = "#4e8ce4,#a6c6f2,#a6a6a6,#8c8c8c",
                    dest    = "palette",
                    metavar = "",
                    help    = "colour palette [format: <col1,col2,...>]")
parser$add_argument("-r", "--response-lab",
                    type    = "character",
                    dest    = "response_label", 
                    default = "Response",
                    metavar = "",
                    help    = "response label [default: 'Response']")
args <- parser$parse_args()

# Functions -------------------------------------------------------------------

# Function for fixing data
fix_data <- function(data) {

    # Find treatment and response columns and rename then
    names(data)[grep(args$treatment, names(data))] <- "treatment"
    names(data)[grep(args$response, names(data))] <- "response"

    # Remove rows where treatment is missing
    data <- data[data$treatment != "", ]

    # Make into factors
    data$treatment <- factor(data$treatment,
                             levels = sort(unique(data$treatment)))

    # Return fixed data
    return(data)
}

# Function for getting treatment groups
get_groups <- function(data) {

    # Coerce into data frame an get treatment groups
    tukey <- as.data.frame(tukey_test$treatment)
    tukey$treatment <- row.names(tukey)
    names(tukey) <- c("diff", "lwr", "upr", "p.adj", "treatment")

    # Separate treatment groups into two columns
    tukey <- separate(data = tukey,
                      sep  = "-",
                      col  = treatment,
                      into = c("treat_1", "treat_2"))
    tukey <- arrange(tukey, desc(treat_2), desc(treat_1))
    tukey$treatment <- paste(tukey$treat_1, tukey$treat_2, sep = " - ")
    tukey$treatment <- factor(tukey$treatment,
                              levels = unique(tukey$treatment))

    # Add colour groups based on significance
    tukey$colour_groups <- "A"
    tukey[tukey$p.adj <= 1 - args$conf_level, "colour_groups"] <- "A"
    tukey[tukey$p.adj > 1 - args$conf_level, "colour_groups"] <- "B"

    # Return Tukey's HSD data
    return(tukey)
}

# Function for finding groups with significant differences
find_differences <- function(data, tukey) {

    # Proportional width of boxplots
    box_width <- length(unique(data$treatment)) * 0.09

    # Find groups with significant differences
    levels <- tukey_test[["treatment"]][, 4]
    labels <- data.frame(multcompLetters(levels)["Letters"])
    labels$treatment <- row.names(labels)
    labels <- arrange(labels, treatment)

    # Separate groups where applicable
    labels <- separate(labels,
                       sep  = 1,
                       col  = "Letters",
                       into = c("group_1", "group_2"))
    labels[labels$group_2 == "", "group_2"] <-
        labels[labels$group_2 == "", "group_1"]

    # Calculate quantiles for fill colours of boxplots
    quantiles <- data %>%
        group_by(treatment) %>%
        summarise(y1 = quantile(response, probs = 0.25),
                  y2 = quantile(response, probs = 0.75))

    # Merge with labels
    shading <- merge(labels, quantiles, by = "treatment")

    # Find coordinates
    shading$x1 <- as.numeric(row.names(shading)) - box_width / 2
    shading$x2 <- as.numeric(row.names(shading)) + box_width / 2

    # Merge with data
    data <- merge(data, shading[c(1,2)], by = "treatment")

    # Initialise empty triangle data objects
    group_2 <- c()
    triangle <- c()
    t2_x =  c()
    t2_y <- c()

    # Get triangle data
    for (n in c(1:nrow(shading))) {
        group_2 <- c(group_2, rep(shading[n, "group_2"], 3))
        triangle <- c(triangle, n, n, n)
        t2_x <- c(t2_x, shading[n, "x1"], shading[n, "x2"], shading[n, "x2"])
        t2_y <- c(t2_y, shading[n, "y1"], shading[n, "y1"], shading[n, "y2"])
    }
    triangle_2 <- data.frame(group    = group_2,
                             triangle = triangle,
                             x        = t2_x,
                             y        = t2_y)

    # Return triangle
    return(list(data, triangle_2))
}

# Function for plotting results
plot_results <- function(data, triangle) {
    
    # Proportional width of boxplots
    box_width <- length(unique(data$treatment)) * 0.09

    # Change underscores to spaces in treatment names
    data$treatment <- gsub("_", " ", data$treatment)

    # Get colour palette
    palette <- strsplit(args$palette, ",")[[1]]

    # Boxplot
    gg_box <- ggplot(data, aes(x    = treatment,
                               y    = response,
                               fill = group_1)) +
        geom_boxplot(outlier.shape = NA,
                     colour        = "#4d4d4d",
                     width         = box_width) +
        geom_polygon(data = triangle, aes(x     = x,
                                          y     = y,
                                          group = triangle,
                                          fill  = group)) +
        geom_boxplot(outlier.shape = NA,
                     colour        = "#4d4d4d",
                     width         = box_width, aes(fill = NA)) +
        stat_boxplot(data   = data,
                     geom   = "errorbar",
                     width  = 0.15, 
                     colour = "#4d4d4d",
                     aes(x = treatment,
                         y = response)) +
        geom_point(position = position_jitter(w = 0.15, h = 0),
                   stroke   = 0.5,
                   alpha    = 0.75, 
                   size     = 1,
                   colour   = "#4d4d4d") +
        theme_bw() +
        theme(legend.position = "none") +
        labs(x = NULL,
             y = args$response_label) +
        scale_fill_manual(values = palette)

    # Confidence interval plot
    limit <- max(abs(tukey$lwr), tukey$upr) * 1.05
    gg_conf <- ggplot(tukey, aes(x      = treatment,
                                 y      = diff,
                                 ymin   = lwr,
                                 ymax   = upr, 
                                 colour = colour_groups)) +
        coord_flip(ylim = c(-limit, limit)) +
        geom_pointrange(shape = 20,
                        size  = 0.3) +
        geom_errorbar(width = 0.3) +
        theme_bw() +
        geom_hline(yintercept = 0,
                   colour = "#4d4d4d",
                   linetype = 2) +
        theme(panel.grid.major.x = element_blank(), 
              panel.grid.minor.x = element_blank(),
              axis.text          = element_text(size = 7),
              plot.title         = element_text(hjust = 0.5),
              legend.position    = "none") +
        labs(y = NULL , x = NULL) +
        scale_colour_manual(values = c(palette[1], "#a6a6a6"))

    # Arrange in grid
    gg <- grid.arrange(gg_box, gg_conf, ncol = 2)

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
suppressPackageStartupMessages(library("multcompView"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))

# Fix data for ANOVA
data <- fix_data(data)

# Perform ANOVA
message("Performing ANOVA and Tukey's HDS ...")
model <- lm(response ~ treatment, data = data)
aov <- aov(model)

# Perform Tukey's Honest Significant Difference test
tukey_test <- TukeyHSD(aov, conf_level = args$conf_level)

# Get treatment groups
tukey <- get_groups(tukey_test)

# Check for significant differences and prepare for plotting
plot_prep <- find_differences(data, tukey)
data <- plot_prep[[1]]
triangle <- plot_prep[[2]]

# Plot
message("Plotting results ...")
gg <- plot_results(data, triangle)

# Save to file
ggsave(args$output, gg, dpi = 300, width = 14, height = 7)
message("Done.")
