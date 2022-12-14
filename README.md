# INTACT

INTACT (INtegration of TWAS And ColocalizaTion) is an R package for constraining 
transcriptome-wide association study correlation statistics through an empirical 
Bayes approach. As input, INTACT takes gene-level colocalization probabilities 
and TWAS z-statistics, returning a posterior probability of putative causality 
for each gene.

The package also performs gene set enrichment estimation using probabilistic 
INTACT results and a list of pre-defined gene sets (INTACT-GSE).

For a thorough description of the methodology, refer to
https://doi.org/10.1101/2022.07.19.500651


## Installing via Github

```
library(devtools)
devtools::install_github("jokamoto97/INTACT")
```

## Installing via Bioconductor

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("INTACT")
```

## Examples

Here is an example of how to integrate colocalization probabilites and TWAS 
z-statistics for a simulated data set ```simdat```:

```
intact(GLCP_vec=simdat$GLCP, z_vec = simdat$TWAS_z)
```

To perform gene set enrichment analysis using simdat and pre-defined gene set 
list ```gene_set_list```, run

```
intactGSE(gene_data = simdat, gene_sets = gene_set_list)
```


## Support

Please contact xwen@umich.edu or jokamoto@umich.edu if you have any questions.
