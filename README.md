# Miscellaneous Bioinformatic Tools

This is a collection of miscellaneous bioinformatic tools and scripts that I
have created and used during my PhD studies. They are written in R, Python and
bash (but mostly R) and are meant to be called from the command line. They are
sorted according to their general theme and use. You can get help by calling
any script from the command line with the `-h` (or `--help`) flag.

## Variant analysis with *seqCAT*

These are scripts using the [seqCAT][seqcat] Bioconductor R-package for variant
analysis of high throughput sequencing (HTS) data. They are simple wrappers for
a common workflow for seqCAT, to facilitate re-use and aggregate analysis of
many samples at once.

```{bash Variant analysis}
# Create SNV profiles from VCFs in a directory
create_profiles.R <input VCF directory> <output profile directory>

# Compare genetic similarities for all created profiles
compare_profiles.R <input profile directory> similarities.txt

# Create a similarity heatmap from comparisons
similarity_heatmap.R similarities.txt similarities.png --cluster
```

## Get gene and transcript information from biomaRt

This is a simple script for fetching information regarding gene and transcript
IDs, symbols and biotypes from Ensembl using the [biomaRt][biomart] package.
These information files are used in some of the other scripts for filtering
out non-coding genes or mapping transcripts to genes.

```{bash Get biomaRt info}
# Get biomaRt information
get_biomart_info.R GRCh38 biomart/biomart.grch38.txt
```

## Differential expression analysis

These are scripts related to differential expression analyses (DEA) of RNA-seq
data. The main one, `de_analysis.R` calculates which genes are differentially
expressed between two samples, using either [DESeq2][deseq2], [edgeR][edger] or
[Limma][limma] or a combination of overlapping DEGs between them. It can handle
inputs as either raw counts or gene expression measurements (TPM) from both
[Kallisto][kallisto] and [Salmon][salmon], using [TXimport][tximport]. The
other scripts can analyse the output from `de_analysis.R` in various ways.

```{bash DEA}
# Calculate differential expression
# This script has many options and parameters: use `-h` to see them all
de_analysis.R counts.txt sample1,sample2 biomaRt_info.txt degs.txt edgeR

# Plot the p-value distribution of DEGs
p_distribution.R degs.txt p_distribution.png

# Create a volcano plot of DEGs
volcano.R degs.txt volcano.png
```

The `pathway_analysis.R` script can analyse a DEG list and one or more
specified KEGG pathways, finding which genes are differentially expressed
in that pathway and visualises them in an image with fold change colour
gradients. It can also quantify and list the various types of interactions
(*e.g.* phosphorylations, direct interactions, etc.) in said pathway(s),
where a *perturbation event* is defined as an interaction **A --> B** where
**A** is a DEG. All this is done via the [Pathview][pathview] package.

```{bash Pathway analysis}
# Analyse the MAPK pathway
pathway_analysis.R degs.txt
```

[biomart]: https://bioconductor.org/packages/release/bioc/html/biomaRt.html
[deseq2]: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[edger]: http://bioconductor.org/packages/release/bioc/html/edgeR.html
[kallisto]: https://pachterlab.github.io/kallisto/
[limma]: http://bioconductor.org/packages/release/bioc/html/limma.html
[pathview]: https://bioconductor.org/packages/release/bioc/html/pathview.html
[salmon]: https://combine-lab.github.io/salmon/
[seqcat]: https://bioconductor.org/packages/release/bioc/html/seqCAT.html
[tximport]: https://bioconductor.org/packages/release/bioc/html/tximport.html
