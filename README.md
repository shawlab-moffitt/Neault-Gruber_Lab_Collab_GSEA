# Neault-Gruber_Lab_Collab_GSEA

## GSEA Overview

Gene Set Enrichment Analysis (GSEA) is a computational method which examines a set or multiple sets of genes to determine statistical significance differences between two biological states. These gene sets can be defined by various signature databases, such as MSigDB or can be self defined depending on genes of interest.

**Input files:** when running GSEA in R, one should include a \*.gmt (Gene Matrix Transposed) file and a file containing gene expression data for each sample. The GMT file will contain one gene set per row, columns including the gene set name, an optional description column or filled in 'na', followed by a list of genes (can be unequal). The gene expression data should contain 2 groups with differing biological states with the first column or row names being gene names, followed by columns named by sample containing the expression value for each gene.

## ssGSEA Overview

As an extension of GSEA, Single Sample Gene Set Enrichment Analysis (ssGSEA), calculates the enrichment scores seperately for each samples and gene set. The results represent if the particular gene set is overall up- or down-regulated for an individual sample.

## Our Analysis

The expression data for each gene of the phenotypes used in this analysis were ranked by signal-to-noise ratio. The signal-to-noise metric requires at least two catagorical phenotypes with at least three samples for each phenotype. This ratio uses the difference of means by the standard deviation, so the larger the ratio, the larger the difference of means, representing the more distinct the gene expression is in each phenotype.

The GSEA() function from the ClusterProfiler package is used to integrate the ranked gene expression with the gene sets to determine the enrichment scores and P-values. Enrichment plots can then be generated based off of individual gene sets through the gseaplot2() function from the enrichplot package.



