#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "pathview", "XML")
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

# Command parser
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(epilog = paste("Visualise DEGs in chosen pathways,",
    "with gene names (ENSEMBL IDs) and fold change values. Output pathway is",
    "by default the KEGG MAPK/ERK pathway, but several can be specified.", 
    "Output is by default single-layered gene names as per KEGG standard."))
parser$add_argument("input",
                    type    = "character", 
                    help    = "input file path(s) [format: <i1,i2,...,iN>]")
parser$add_argument("-l", "--layers",
                    action  = "store_false",
                    dest    = "layers", 
                    help    = "two-layered image (official gene symbols)")
parser$add_argument("-i", "--image_only",
                    action  = "store_true", 
                    dest    = "image_only",
                    help    = "no XML analysis, only image")
parser$add_argument("-G", "--gene-list",
                    type    = "character", 
                    dest    = "gene_list",
                    default = "",
                    metavar = "",
                    help    = "save pathway genes as list [biomart info file]")
parser$add_argument("-p", "--pathway",
                    type    = "character",
                    dest    = "pathway",
                    default = "04010", 
                    metavar = "",
                    help    = "KEGG IDs [default: 04010; format: <p1,p2,...>]")
parser$add_argument("-P", "--prefix",
                    type    = "character",
                    dest    = "prefix",
                    default = "",
                    metavar = "",
                    help    = "prefix to output file(s)")
parser$add_argument("-c", "--colour",
                    type    = "character",
                    default = "#e4bf4e,#4e8ce4",
                    dest    = "colours",
                    metavar = "",
                    help    = "fold change colours [format: <high,low>]")
parser$add_argument("-r", "--range",
                    type    = "integer",
                    dest    = "range",
                    default = 5, 
                    metavar = "",
                    help    = "range for log2FC colours (default: 5)")
parser$add_argument("-b", "--bins",
                    type    = "integer",
                    dest    = "bins",
                    default = 20, 
                    metavar = "",
                    help    = "number of colour bins (default: 20)")
parser$add_argument("-s", "--sum-method",
                    type    = "character",
                    dest    = "sum_method",
                    default = "mean",
                    metavar = "",
                    help = "sum, mean, median, max; [default: mean]")
args <- parser$parse_args()

# Functions -------------------------------------------------------------------

# Function for writing pathway gene lists to file
write_gene_list <- function(genes, biomart_info) {

    # Get list of genes in pathway
    gene_list <- sort(unlist(strsplit(as.character(genes$name),", ")))
    gene_list <- data.frame(gene = gsub("\\.\\.\\.", "", gene_list))

    # Read biomaRt info
    info <- read.table(args$biomart,
                       sep              = "\t",
                       header           = TRUE,
                       stringsAsFactors = FALSE)
    info <- unique(info[c("ensembl_gene_id", "hgnc_symbol")])

    # Add ENSGIDs and write to file
    gene_list <- merge(gene_list,
                       info,
                       by.x = "gene",
                       by.y = "hgnc_symbol")
    write.table(gene_list, paste0(pathway_name, ".gene_list.txt"), 
                sep       = "\t",
                row.names = FALSE)
}

# Function for analysing pathways XMLs
analyse_pathway_xml <- function(doc) {

    # Get genes and labels
    genes <- XML:::xmlAttrsToDataFrame(
        doc["/pathway/entry[@type = 'gene']/graphics"])
    labels <- XML:::xmlAttrsToDataFrame(
        doc["/pathway/entry[@type = 'gene']"])

    # Write pathway genes to file (if applicable)
    if (args$gene_list != "") {
        write_gene_list(genes, args$gene_list)
    }

    # Get graphic names and merge ID/label dataframes
    genes$graphic.name <- sapply(strsplit(as.character(genes$name),","),
                                 "[", 1)
    mapping <- merge(genes, labels, by = "row.names", sort = FALSE)
    mapping <- mapping[c("id", "graphic.name")]

    # Create empty data frame and loop over relation nodes
    interactions <- data.frame()
    for (node in xpathApply(doc, "/pathway/relation")) {
      
        # Get entries and subtypes
        entries <- t(data.frame(xmlAttrs(node), stringsAsFactors = FALSE))
        subtypes <- data.frame(sapply(xmlChildren(node), xmlAttrs), 
                    stringsAsFactors = FALSE)

        # Merge
        node.m <- merge(entries, subtypes, all = TRUE)[1, ]

        # Collect data in "interactions" dataframe
        interactions <- merge(interactions, node.m, all = TRUE)
    }

    interactions$subtype[!is.na(interactions$subtype.1)] <- 
        paste(interactions$subtype.1[!is.na(interactions$subtype.1)],
              interactions$subtype[!is.na(interactions$subtype.1)], sep = ": ")
    interactions <- interactions[c(1:4)]

    # Initiate pathways stats with total number of pathway interactions
    path_stats <- aggregate(interactions$subtype,
                            by  = list(interactions$subtype),
                            FUN = function(x) sum(!is.na(x)))
    names(path_stats) <- c("Interaction", "Total pathway interactions")

    # Convert interaction entry IDs to graphic names
    interactions <- merge(interactions, mapping,
                          by.x = "entry1",
                          by.y = "id", 
                          sort = FALSE)
    interactions <- merge(interactions, mapping,
                          by.x = "entry2",
                          by.y = "id", 
                          sort = FALSE)

    # Get fold change for in/out interactions
    ge1_unique <- unique(pv_out$plot.data.gene[c("labels", "ge1")])
    interactions <- merge(interactions, ge1_unique,
                          by.x = "graphic.name.x", 
                          by.y = "labels",
                          all.x = TRUE,
                          sort = FALSE)
    interactions <- merge(interactions, ge1_unique,
                          by.x = "graphic.name.y", 
                          by.y = "labels",
                          all.x = TRUE,
                          sort = FALSE)
    interactions <- interactions[c("graphic.name.x", "ge1.x", "subtype")]
    names(interactions) <- c("Gene", "log2FC", "Interaction")

    # Initiate pathway statistics dataframe
    events <- aggregate(interactions$log2FC,
                        by  = list(interactions$Interaction), 
                        FUN = function(x) sum(!is.na(x)))
    names(events) <- c("Interaction" , "total events")
    path_stats <- merge(path_stats, events,
                        by  = "Interaction",
                        all = TRUE)

    # Count up-regulated genes
    if (nrow(subset(interactions, log2FC > 0, select = log2FC)) > 0) {
        events_up <- aggregate(subset(interactions, log2FC > 0,
                                      select = log2FC), 
                                      by = subset(interactions, log2FC > 0,
                                                  select = Interaction),
                      FUN = function(x) sum(!is.na(x)))
    } else {
        events_up <- data.frame("interaction" = NA, "events_up" = NA)
    }

    # Count down-regulated genes
    if (nrow(subset(interactions, log2FC < 0, select = log2FC)) > 0) {
        events_down <- aggregate(subset(interactions, log2FC < 0,
                                        select = log2FC), 
                                        by = subset(interactions, log2FC < 0,
                                                    select = Interaction), 
                        FUN = function(x) sum(!is.na(x)))
    } else {
        events_down <- data.frame("interaction" = NA, "events_down" = NA)
    }

    # Merge with stats
    regulation <- merge(events_up, events_down,
                        by  = "Interaction",
                        all = TRUE)
    path_stats <- merge(path_stats, regulation,
                        by  = "Interaction",
                        all = TRUE)

    # Set final column names
    names(path_stats) <- c("Interaction",
                           "Total pathway interactions", 
                           "Total events",
                           "Up-regulated events",
                           "Down-regulated events")

    # Total column values, naming, etc.
    path_stats <- rbind(path_stats, c(NA, round(colSums(path_stats[,-1], 
                                                        na.rm = TRUE), 0)))
    path_stats[path_stats == 0] <- NA
    path_stats <- path_stats[apply(path_stats, 1,
                                   function(x) any(!is.na(x))), ]
    path_stats$Interaction[length(path_stats$Interaction)] <- "Total"
    path_stats$Interaction <- paste0(toupper(
        substring(path_stats$Interaction, 1, 1)),
        substring(path_stats$Interaction, 2))

    # Write interactions and stats to separate files
    write.table(interactions,
                paste0(args$prefix, pathway_name, ".list.txt"),
                row.names = FALSE,
                sep = "\t",
                na = "")
    write.table(path_stats,
                paste0(args$prefix, pathway_name, ".stats.txt"),
                row.names = FALSE,
                sep = "\t",
                na = "")
}

# Analysis --------------------------------------------------------------------

# Read data
message("Reading data ...")
inputs <- strsplit(args$input, ",")[[1]]
data <- read.table(inputs[1],
                   sep              = "\t",
                   row.names        = 1,
                   header           = TRUE,
                   fill             = TRUE,
                   stringsAsFactors = FALSE)

# Read all other input files (if applicable)
if (length(inputs) > 1) {
    for (n in 2:length(inputs)) {

        # Read data
        temp <- read.table(inputs[n],
                           sep              = "\t",
                           row.names        = 1,
                           header           = TRUE,
                           fill             = TRUE,
                           stringsAsFactors = FALSE)

        # Merge with existing data
        data <- merge(data, temp, by = "row.names", all = TRUE)
        row.names(data) <- data$Row.names
        data$Row.names <- NULL
    }
}

# Load pathway analysis packages
message("Performing pathway analysis ...")
suppressPackageStartupMessages(suppressWarnings(library("pathview")))
suppressPackageStartupMessages(suppressWarnings(library("XML")))

# Get pathways to analyse
pathways <- strsplit(args$pathway, ",")[[1]]

# Get colours
colours <- strsplit(args$colours, ",")[[1]]

# id_mapping
id_map_entrez <- suppressMessages(id2eg(ids      = row.names(data),
                                        category = "ENSEMBL",
                                        org      = "Hs"))
data_entrez <- mol.sum(id.map   = id_map_entrez,
                       mol.data = data[grep("log2FC",
                                            names(data),
                                            value = TRUE)])

# Loop through each pathway
for (pathway in pathways) {
  
  # Pathview
    pv_out <- suppressMessages(pathview(gene.data   = data_entrez,
                                        pathway.id  = pathway,
                                        species     = "hsa",
                                        kegg.native = TRUE,
                                        map.symbol  = TRUE,
                                        same.layer  = args$layers,
                                        node.sum    = args$sum_method,
                                        limit       = args$range,
                                        bins        = args$bins,
                                        low         = c(rep(colours[2], 2)), 
                                        mid         = c(rep("white", 2)),
                                        high        = c(rep(colours[1], 2))))
    # Read pathway XML
    doc <- xmlParse(paste("hsa", pathway, ".xml", sep = ""), 
                    useInternalNodes = TRUE)

    # pathway_name, gene names and entry IDs
    pathway_name <- strsplit(xpathApply(doc, "/pathway", xmlGetAttr, 
                                        "title")[[1]], " ")[[1]][1]

    # XML analysis for single input (if applicable)
    if (all(!args$image_only, !(length(inputs) > 1)) ) {
        message(paste("Analysing", pathway_name, "pathway ..."))
        analyse_pathway_xml(doc)
    }

    # Remove superfluous pathway files
    file.remove(paste("hsa", pathway, ".xml", sep = ""))
    file.remove(paste("hsa", pathway, ".png", sep = ""))

    # Rename pathway image
    name <- Sys.glob(paste("*", pathway, "*\\.png", sep = ""))
    file.rename(name, paste0(args$prefix, pathway_name, ".png"))
}

# End message
message("Done.")
