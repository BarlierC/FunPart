![funpart package](pictures/FunPart_logo.png) 

# FunPart

![funpart version](https://img.shields.io/static/v1?label=funpart&message=v1.0&color=green)

FunPart: Deciphering and characterizing functional heterogeneity in single cell data

FunPart is a computational tool that partitions heterogeneous cell populations into functionally distinct subpopulations and simultaneously identifies modules of functionally relevant set of genes for each of them.

## INSTALLATIONS

#### Installation of devtools from CRAN

install.packages("devtools")

library("devtools")

#### Installation of FunPart R package using Devtools




### EXAMPLES

```R

# Load the packages and dependencies

require(funpart)
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

#Set seed
set.seed(1234)
clustnum <- parallel::makeCluster(5)
doParallel::registerDoParallel(clustnum)

options(stringsAsFactors = FALSE)


# Load data necessary to run FunPart
# Mouse Biological Processes (GO) annotation file (gda)
data(gda)
# Mouse TFs 
data(mouse_tfs)

#Load your data: cells in columns, genes in rows
#data_exp

#Run the Funtional Splitting with default parameters
res <- run_functional_splitting(data_exp,mouse_tfs$Symbol,gda)

```
