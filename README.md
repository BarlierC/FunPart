![funpart package](pictures/FunPart_logo.png) 

# FunPart

![funpart version](https://img.shields.io/static/v1?label=funpart&message=v1.0&color=green) ![funpart licence](https://img.shields.io/badge/licence-GPL-blue)

#### FunPart: Deciphering and characterizing functional heterogeneity in single cell data

FunPart is a computational tool that partitions heterogeneous cell populations into functionally distinct subpopulations and simultaneously identifies modules of functionally relevant set of genes for each of them.

## HOW TO INSTALL IT

#### Installation of FunPart using Devtools

```R
install.packages("devtools")
library("devtools")
install_github(“BarlierC/FunPart”)
```

## HOW TO USE IT

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
