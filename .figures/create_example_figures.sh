#!/bin/bash

# Differential expression figures
p_distribution.R data/degs.txt example_p_distribution.png
volcano.R data/degs.txt example_volcano.png
enrichment_analysis.R data/degs.txt \
    ../biomart/biomart.grch38.txt \
    example_enrichment.png -t GOslim
pathway_analysis.R data/degs.txt
mv MAPK.png example_pathway.png
rm MAPK*

# Expression figures
expression_barplot.R data/tpm_transcripts.txt ENST00000275493,ENST00000420316 \
    example_transcripts.png -s sample_1,sample_2

# Unsupervised learning figures
pc_analysis.R data/tpm_genes.txt example_pca.png -l -t
pairwise_correlations.R data/tpm_genes.txt data/temp.txt -l
add_metadata.R data/temp.txt data/metadata.txt data/corr_genes.txt
rm data/temp.txt
heatmap.R data/corr_genes.txt r2 example_heatmap.png -c all -g group
dendrogram.R data/corr_genes.txt r2 example_dendrogram.png -k 5 -g group
mds.R data/corr_genes.txt r2 example_mds.png -g group

# Misc figures
anova.R data/anova_data.txt example_anova.png treatment response
