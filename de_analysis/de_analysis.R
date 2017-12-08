#!/usr/bin/env Rscript

# Install missing packages (if applicable)
packages <- c("argparse", "tximport", "rhdf5", "DESeq2", "edgeR", "limma")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  message("installing missing packages ...")
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

# Command parser --------------------------------------------------------------

suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser(epilog = paste("Calculates differential",
    "expression using various methods (DESeq2, edgeR and Limma) or their",
    "overlaps. Input data can be either raw counts or TPM from Salmon",
    "or Kallisto."))
parser$add_argument("input",
                    type    = "character",
                    help    = "input file path")
parser$add_argument("samples",
                    type    = "character",
                    help    = "samples to analyse [format: <case,control>]")
parser$add_argument("biomart",
                    type    = "character",
                    help    = "path to biomaRt info file for genes / biotypes")
parser$add_argument("output",
                    type    = "character",
                    help    = "output file path")
parser$add_argument("mode",
                    type    = "character", 
                    help    = paste("[intersection / partial / union /",
                                    "DESeq / edgeR / Limma"))
parser$add_argument("-n", "--non-coding",
                    action  = "store_true",
                    dest    = "non_coding",
                    help    = "also include non-coding genes in the analysis")
parser$add_argument("-N", "--non_degs",
                    action  = "store_true",
                    dest    = "non_degs",
                    help    = "include non-DEGs in results [default: FALSE]")
parser$add_argument("-f", "--FDR-cutoff",
                    type    = "double",
                    dest    = "fdr_cutoff",
                    metavar = "",
                    default = 0.01,
                    help    = "desired FDR cutoff [default: 0.01]")
parser$add_argument("-F", "--fc-cutoff",
                    type    = "double",
                    dest    = "fc_cutoff",
                    metavar = "",
                    default = 2,
                    help    = "desired fold change cutoff [default: 2]")
parser$add_argument("-r", "--raw-cutoff",
                    type    = "integer",
                    dest    = "raw_cutoff", 
                    metavar = "",
                    default = 0,
                    help    = "cutoff for input data")
parser$add_argument("-e", "--expression",
                    type    = "character",
                    dest    = "expression", 
                    metavar = "",
                    default = "",
                    help    = "input is a [kallisto/salmon] directory")
parser$add_argument("-i", "--input-2",
                    type    = "character",
                    dest    = "input_2",
                    metavar = "",
                    default = "",
                    help    = "second input file path")
args = parser$parse_args()

# Function definitions --------------------------------------------------------

# Function for reading biomaRt info
read_biomart_info <- function(biomart_file, non_coding) {

    # Read info file
    message("Reading biomaRt info ...")
    info <- read.table(biomart_file,
                               sep              = "\t",
                               header           = TRUE,
                               stringsAsFactors = FALSE)

    # Remove non-coding genes (if applicable)
    biomart_info <- info[info$gene_biotype == "protein_coding", ]

    # Return biomaRt info dataframe
    return(biomart_info)
}

# Function for reading count data
read_counts <- function(file_1, file_2, sample_1, sample_2) {

    # Read data
    message("Reading count data ...")
    data <- read.table(file_1, header = TRUE, row.names = 1)

    if (file_2 == "") {

        # Read both samples from same input file
        data_1 <- data[, c(grep(sample_1, names(data)))]
        data_2 <- data[, c(grep(sample_2, names(data)))]
        names(data_1) <- paste("x.", 1:length(names(data_1)), sep = "")
        names(data_2) <- paste("y.", 1:length(names(data_2)), sep = "")

    } else {

        # Read first input file for sample 1
        data_1 <- data[, c(grep(samples[1], names(data)))]
        names(data_1) <- paste("x.", 1:length(names(data_1)), sep = "")

        # Read second input file for sample 2
        data_2 <- read.table(file_2, header = TRUE, row.names = 1)
        data_2 <- data_2[, c(grep(samples[2], names(data_2)))]
        names(data_2) <- paste("y.", 1:length(names(data_2)), sep = "")
    }

    # Merge data
    data <- merge(data_1, data_2, by = "row.names")
    row.names(data) <- data$Row.names
    data$Row.names <- NULL

    # Return data
    return(data)
}
 
# Function for reading Kallisto / Salmon expression data
read_expression <- function(dir, software, biomart_info, sample_1, sample_2) {

    # Load packages
    message("Reading expression data ...")
    suppressPackageStartupMessages(suppressWarnings(library("tximport")))

    # Define file types
    if (tolower(software) == "kallisto") {
        file_type <- "abundance.h5"
    } else if (tolower(software) == "salmon") {
        file_type <- "quant.sf"
    }

    # Find all appropriate files for each sample
    files <- list.files(dir,
                        pattern    = file_type,
                        recursive  = TRUE,
                        full.names = TRUE)
    files_1 <- grep(sample_1, files, value = TRUE)
    files_2 <- grep(sample_2, files, value = TRUE)
    paths <- c(files_1, files_2)
    names(paths) <- c(paste("x", 1:sum(grepl(sample_1, paths)), sep = "."), 
                      paste("y", 1:sum(grepl(sample_2, paths)), sep = "."))

    # Get transcript and gene info
    tx2gene <- biomart_info[c("ensembl_transcript_id", "ensembl_gene_id")]

    # Import abundance estimates
    data <- suppressMessages(tximport(paths,
                                      type                = software,
                                      tx2gene             = tx2gene,
                                      countsFromAbundance = "lengthScaledTPM"))
    data <- as.data.frame(data$counts)

    # Return data
    return(data)
}

# Function for defining DE conditions
define_conditions <- function(data) {

    # Define the DE conditions
    n_1 <- length(grep("x.", names(data)))
    n_2 <- length(grep("y.", names(data)))
    condition <- as.factor(c(rep("x", n_1),
                             rep("y", n_2)))
    condition <- factor(condition, levels = c("y", "x"))

    # Return conditions
    return(condition)
}

# Function for finalising DAE: add gene names, regulation and non-degs
finalise_dea <- function(res, biomart_info, non_degs, fdr_cutoff, fc_cutoff) {
    
    # Get gene-only ID list
    genes <- biomart_info[c("ensembl_gene_id", "hgnc_symbol")]
    genes <- genes[!duplicated(genes$ensembl_gene_id), ]
    names(genes) <- c("ENSGID", "gene")

    # Merge gene names with DEA results
    res <- merge(res,
                 genes,
                 by.x = "row.names", 
                 by.y = "ENSGID",
                 all.x = TRUE)
    row.names(res) <- res$Row.names
    res$Row.names <- NULL

    # Remove non-degs (if applicable)
    if(!non_degs) {
        res <- res[res$FDR <= fdr_cutoff &
                   !is.na(res$FDR) &
                   2^abs(res$log2FC) >= fc_cutoff, ]
    }

    # Check and add regulation direction
    res[res$log2FC >= 1, "direction"] <- "Up"
    res[res$log2FC <= -1, "direction"] <- "Down"

    # Return data
    return(res)
}

# Function for DESeq2 analysis
dea_deseq <- function(data, biomart_info, non_degs, fdr_cutoff, fc_cutoff) {

    # Load package
    message("Performing DESeq2 analysis ...")
    suppressPackageStartupMessages(suppressWarnings(library("DESeq2")))

    # Get DEA conditions
    condition <- define_conditions(data)

    # DE calculations
    organization <- data.frame(row.names = names(data), condition = condition)
    dds <- suppressMessages(DESeqDataSetFromMatrix(countData = round(data, 0),
                                                   colData   = organization, 
                                                   design    = ~condition))
    res <- data.frame(results(suppressMessages(DESeq(dds))))

    # Change column names and finalise DEA
    names(res) <- c("baseMean", "log2FC", "lfcSE", "stat", "p_value", "FDR")
    res <- finalise_dea(res, biomart_info, non_degs, fdr_cutoff, fc_cutoff)

    # Return DESeq2 results
    return(res)
}

# Function for edgeR DEA
dea_edger <- function(data, biomart_info, non_degs, fdr_cutoff, fc_cutoff) {

    # Load packages
    message("Performing edgeR analysis ...")
    suppressPackageStartupMessages(suppressWarnings(library("edgeR")))

    # Define DEA conditions
    condition <- define_conditions(data)

    # Estimate disperions
    dge_edger <- DGEList(counts = data, group = condition)
    dispersion <- estimateCommonDisp(dge_edger)
    dispersion <- estimateTagwiseDisp(dispersion)

    # Exact test and  results table
    et <- exactTest(dispersion)
    res <- data.frame(row.names = row.names(et$table), c(et$table))
    names(res) <- c("log2FC", "logCPM", "p_value")

    # Multiple comparison adjustment and significance testing
    res$FDR <- p.adjust(res$p_value, method = "BH")
    
    # Finalise DEA
    res <- finalise_dea(res, biomart_info, non_degs, fdr_cutoff, fc_cutoff)

    # Return edgeR results
    return(res)
}

# Function for Limma DEA
dea_limma <- function(data, biomart_info, non_degs, fdr_cutoff, fc_cutoff) {

    # Load packages
    message("Performing Limma analysis ...")
    suppressPackageStartupMessages(suppressWarnings(library("limma")))
    suppressPackageStartupMessages(suppressWarnings(library("edgeR")))

    # Get DEA conditions
    condition <- define_conditions(data)

    # TMM normalization and voom transformation
    dge_limma <- calcNormFactors(DGEList(counts = data))
    design <- model.matrix(~factor(condition))
    v <- voom(dge_limma, design, plot = FALSE)

    # Perform DEA
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    res <- suppressMessages(topTable(fit, number = nrow(fit)))
    res <- res[order(row.names(res)), ]

    # Change column names and finalise DEA
    names(res) <- c("log2FC", "aveExpr", "t", "p_value", "FDR", "B")
    res <- finalise_dea(res, biomart_info, non_degs, fdr_cutoff, fc_cutoff)

    # Return Limma results
    return(res)
}

# Function for calculating DEA overlaps
dea_overlaps <- function(res_1, res_2, res_3, mode) {

    # Calculate the union of all three sets
    message("Calculating gene overlaps ...")
    res_union <- merge(res_1[c("gene", "log2FC")],
                       res_2[c("gene", "log2FC")],
                       by  = "row.names",
                       all = TRUE)
    res_union <- merge(res_union,
                       res_3[c("gene", "log2FC")],
                       by.x = "Row.names",
                       by.y = "row.names",
                       all  = TRUE)

    # Calculate union fold change means
    row.names(res_union) <- res_union$Row.names
    res_union$FC <- rowMeans(res_union[, c("log2FC", "log2FC.x", "log2FC.y")],
                             na.rm = TRUE)
    res_union <- data.frame(ENSGID = row.names(res_union),
                            gene   = res_union$gene,
                            log2FC = res_union$FC)

    # Add regulation direction
    res_union[res_union$log2FC >= 1, "direction"] <- "Up"
    res_union[res_union$log2FC <= -1, "direction"] <- "Down"

    # Calculate intersections (if applicable)
    if (mode == "INTERSECTION") {

        # Calculate total intersection
        inter_table <- table(c(row.names(res_1),
                               row.names(res_2), 
                               row.names(res_3)))
        res_inter <- data.frame(ENSGID = names(inter_table)[inter_table == 3])
        res_inter <- merge(res_inter, res_union)

        # Return total intersection
        return(res_inter)

    } else if (mode == "PARTIAL") {
        
        # Calculate partial intersection
        inter_table <- table(c(row.names(res_1),
                               row.names(res_2), 
                               row.names(res_3)))
        res_part <- data.frame(ENSGID = names(inter_table)[inter_table  >= 2])
        res_part <- merge(res_part, res_union)

        # Return partial intersection
        return(res_part)

    } else if (mode == "UNION") {

        # Return union
        return(res_union)
    }
}

## Analysis -------------------------------------------------------------------

# Check for valid mode
mode <- toupper(args$mode)
modes <- c("INTERSECTION", "PARTIAL", "UNION", "DESEQ2", "EDGER", "LIMMA")
if (!mode %in% modes) {
  stop(paste0("invalid run mode, please use one of the following:\n",
              "INTERSECTION, PARTIAL, UNION, DESeq2, edgeR, Limma"))
}

# Get samples
samples <- strsplit(args$samples, ",")[[1]]
sample_1 <- samples[1]
sample_2 <- samples[2]

# Read biomaRt info
biomart_info <- read_biomart_info(args$biomart)

# Read data
if (args$expression != "") {

    # Data is in a Kallisto / Salmon directory
    data <- read_expression(args$input, args$expression, biomart_info,
                            sample_1, sample_2)

} else {

    # Data is in a counts file
    data <- read_counts(args$input, args$input_2, sample_1, sample_2)
}
    
# Remove non-coding genes (if applicable)
if (!args$non_coding) {
    data <- data[row.names(data) %in% biomart_info$ensembl_gene_id, ]
}

# Remove data below cutoff (if applicable)
if (args$raw_cutoff != 0) {

    # Separate samples and remove rows below cutoff
    x <- grep("x.", names(data), value = TRUE)
    y <- grep("y.", names(data), value = TRUE)
    data <- data[rowSums(data[x]) / length(x)  >= args$raw_cutoff |
                 rowSums(data[y]) / length(y)  >= args$raw_cutoff, ]
}

# Run DEA
overlap_modes <- c("INTERSECTION", "PARTIAL", "UNION")
if (mode == "DESEQ2") {  # DESeq2 only

    # Perform DESeq2 DEA
    res <- dea_deseq(data, biomart_info, args$non_degs, args$fdr_cutoff,
                     args$fc_cutoff)

    # Set ENSGID for output
    res$ENSGID <- row.names(res)
    res <- res[c(9,1:8)]

} else if (mode == "EDGER") {  # edgeR only

    # Perform edgeR DEA
    res <- dea_edger(data, biomart_info, args$non_degs, args$fdr_cutoff,
                     args$fc_cutoff)

    # Set ENSGID for output
    res$ENSGID <- row.names(res)
    res <- res[c(7, 1:6)]
    
} else if (mode == "LIMMA") {  # Limma only

    # Perform Limma DEA
    res <- dea_limma(data, biomart_info, args$non_degs, args$fdr_cutoff,
                     args$fc_cutoff)
    
    # Set ENSGID for output
    res$ENSGID <- row.names(res)
    res <- res[c(9, 1:8)]

} else if (mode %in% overlap_modes) {  # Overlap modes

    # Calculate DEA for DESeq2, edgeR and limma
    res_1 <- dea_deseq(data, biomart_info, args$non_degs, args$fdr_cutoff,
                       args$fc_cutoff)
    res_2 <- dea_edger(data, biomart_info, args$non_degs, args$fdr_cutoff,
                       args$fc_cutoff)
    res_3 <- dea_limma(data, biomart_info, args$non_degs, args$fdr_cutoff, 
                       args$fc_cutoff)

    # Calculate overlaps
    res <- dea_overlaps(res_1, res_2, res_3, mode)
}

# Write output to file
write.table(res,
            args$output,
            sep       = "\t",
            row.names = FALSE)
message("Done.")
