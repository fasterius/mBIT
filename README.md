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

[seqcat]: https://bioconductor.org/packages/release/bioc/html/seqCAT.html
