![funpart package](pictures/FunPart_logo.png) 

# FunPart

![funpart version](https://img.shields.io/static/v1?label=funpart&message=v1.0&color=green) ![funpart licence](https://img.shields.io/badge/licence-GPL-blue)

#### FunPart: Deciphering and characterizing functional heterogeneity in single cell data

FunPart is a computational tool that partitions heterogeneous cell populations into functionally distinct subpopulations and simultaneously identifies modules of functionally relevant set of genes for each of them.

## How to install it

#### Installation of FunPart using Devtools

```R
install.packages("devtools")
library("devtools")
install_github("BarlierC/FunPart")
```

## How to use it

```R

# Load packages and dependencies
library(pheatmap)
require(data.tree)
require("plyr")
require("igraph")
require("EnvStats")
require("data.table")
library(stringr)
library(clusterProfiler)
library(Matrix)
library(doParallel)
library(foreach)
require(Seurat)
library(WGCNA)
require(funpart)

#Parallel background, using 4 cores
clustnum <- parallel::makeCluster(4)
doParallel::registerDoParallel(clustnum)

options(stringsAsFactors = FALSE)

# Load necessary data to run FunPart

# 1) The annotation file: annotated immune modules by Singhania et al. (https://doi.org/10.1038/s41467-019-10601-6) 
data(gda)

# 2) Mouse Transcription factors 
data(mouse_tfs)

# 3) Your data (matrix or dataframe with cells in columns and genes in rows)
#data_exp

#Run FunPart with default parameters
res <- run_functional_splitting(data_exp,mouse_tfs$Symbol,gda)

```

## Result object

The output is an **object of class 'functionalSplitting'** composed of **6 slots**:
1. **data** - *contains the filtered and normalized matrix used to perform the analysis*
2. **clust** - *contains the functional cell states identified*
3. **genesets** - *list containing all the genes used to split each level*
4. **cliques** - *list of TFs cliques and genes identified for each branch C1 or C2 and each level*
5. **functionalenrich** - *list of all the significant functional enrichment found for each genes modules*
6. **modules** - *list of all the genes modules identified*
