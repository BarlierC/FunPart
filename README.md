![funpart package](pictures/FunPart_logo.png) 

# FunPart

![funpart version](https://img.shields.io/static/v1?label=funpart&message=v1.0&color=green) ![funpart licence](https://img.shields.io/badge/licence-GPL-blue)

#### FunPart: Deciphering and characterizing functional heterogeneity in single cell data

FunPart is a computational tool that partitions heterogeneous cell populations into functionally distinct subpopulations and simultaneously identifies modules of functionally relevant set of genes for each of them.

## How to install it

#### Installation of FunPart using Devtools

Please note that **FunPart was built under R version 3.5.1**. *A version compatible with R >= 4.0 will be developed*.
*Details about the R environment and packages versions used can be found at the end of the README (R session info).*

```R
install.packages("devtools")
library("devtools")
install_github("BarlierC/FunPart")
```

## How to use it

#### Load all the necessary packages and prepare the parallel environment
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
```

#### Load necessary data to run FunPart
```R
# 1) The annotation file: annotated immune modules by Singhania et al. (https://doi.org/10.1038/s41467-019-10601-6) 
data(gda)

# 2) Mouse Transcription factors 
data(mouse_tfs)

# 3) Your data (matrix or dataframe with cells in columns and genes in rows)
#data_exp
```

#### Run FunPart
```R
#Run FunPart with default parameters
res <- run_functional_splitting(data_exp,mouse_tfs$Symbol,gda)
```

## FunPart parameters
1. **scm** - single cell RNA-seq expression matrix (*cells in columns, genes in rows)
2. **tfs** - vector of TFs names or ID (*it needs to match with the type of genes you have in your matrix*)
3. **gda** - file used to perform the functional enrichment
4. **norm** - perform normalization (*using Seurat, default is TRUE*)
5. **qtarget** - quantile to select the strongest target genes for each identified cliques (*default to 0.90*)
6. **adjMethod** - method used in the enrichment part to perform the p-value adjustement (*for possible values, see enricher package, default = "bonferroni"*)
7. **cutoff** - p-adjusted value cutoff for the functional enrichment (*default=0.05*)
8. **percExp** - percentage of cells a gene needs to be expressed in to be considered (*default=10%)
9. **qExp** - quantile for which the gene expression not be used as considered too low (*default=0.05*)


## FunPart object
The output is an **object of class 'functionalSplitting'** composed of **6 slots**:
1. **data** - *contains the filtered and normalized matrix used to perform the analysis*
2. **clust** - *contains the functional cell states identified*
3. **genesets** - *list containing all the genes used to split each level*
4. **cliques** - *list of TFs cliques and genes identified for each branch C1 or C2 and each level*
5. **functionalenrich** - *list of all the significant functional enrichment found for each genes modules*
6. **modules** - *list of all the genes modules identified*


## Get Summary Table of FunPart object
The function 'getModuleFunctionalState()' will summarize FunPart results into a summary table with 6 columns:
1. **Module**: M1 or M2 (*correspond to cliques slot C1 or C2 of functionalSplitting object respectively*)
2. **Branch**: level and branch (*e.g., 0_0_1 is third level branch 1, 0_1 is second level branch 1 and 0 is first level branch 0*)
3. **Type**: intermediate or direct genes modules. An intermediate genes modules characterize a group of functional states whereas a direct modules genes characterize the specific functional state.
4. **TFs**: transcription factors of the clique (of the genes modules)
5. **Genes**: genes having a direct interaction with TFs of the clique. TFs of the clique + genes having a direct interaction = genes module
6. **Enrichment**: immune processes enriched for the specific module

```R
summary_table <- getModuleFunctionalState(res)
```

## R session info

```R
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X Snow Leopard 11.6

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Seurat_3.1.5           doParallel_1.0.16      iterators_1.0.13       foreach_1.5.1         
 [5] Matrix_1.2-18          clusterProfiler_3.10.1 stringr_1.4.0          data.table_1.14.0     
 [9] EnvStats_2.3.1         igraph_1.2.6           plyr_1.8.6             data.tree_0.7.11      
[13] pheatmap_1.0.12       

loaded via a namespace (and not attached):
  [1] fgsea_1.8.0          Rtsne_0.15           colorspace_2.0-0     ellipsis_0.3.1       ggridges_0.5.2      
  [6] qvalue_2.14.1        leiden_0.3.3         listenv_0.8.0        farver_2.1.0         urltools_1.7.3      
 [11] graphlayouts_0.7.0   ggrepel_0.8.2        bit64_4.0.5          AnnotationDbi_1.44.0 fansi_0.4.2         
 [16] xml2_1.3.2           codetools_0.2-18     splines_3.5.1        cachem_1.0.6         GOSemSim_2.8.0      
 [21] polyclip_1.10-0      jsonlite_1.7.2       ica_1.0-2            cluster_2.1.1        GO.db_3.7.0         
 [26] png_0.1-7            uwot_0.1.8           sctransform_0.2.1    ggforce_0.3.1        BiocManager_1.30.10 
 [31] compiler_3.5.1       httr_1.4.2           rvcheck_0.1.8        lazyeval_0.2.2       assertthat_0.2.1    
 [36] fastmap_1.1.0        tweenr_1.0.1         htmltools_0.5.0      prettyunits_1.1.1    tools_3.5.1         
 [41] rsvd_1.0.3           gtable_0.3.0         glue_1.4.2           RANN_2.6.1           reshape2_1.4.4      
 [46] DO.db_2.9            dplyr_0.8.5          fastmatch_1.1-0      Rcpp_1.0.6           enrichplot_1.2.0    
 [51] Biobase_2.42.0       vctrs_0.3.6          gdata_2.18.0         ape_5.4              nlme_3.1-148        
 [56] lmtest_0.9-38        ggraph_2.0.2         globals_0.14.0       irlba_2.3.3          lifecycle_1.0.0     
 [61] gtools_3.8.2         future_1.18.0        DOSE_3.8.2           zoo_1.8-8            europepmc_0.3       
 [66] MASS_7.3-53.1        scales_1.1.1         tidygraph_1.1.2      hms_0.5.3            RColorBrewer_1.1-2  
 [71] pbapply_1.4-2        reticulate_1.16      memoise_2.0.0        gridExtra_2.3        ggplot2_3.3.3       
 [76] UpSetR_1.4.0         triebeard_0.3.0      stringi_1.4.6        RSQLite_2.2.0        S4Vectors_0.20.1    
 [81] caTools_1.17.1.3     BiocGenerics_0.28.0  BiocParallel_1.16.6  rlang_0.4.10         pkgconfig_2.0.3     
 [86] bitops_1.0-6         lattice_0.20-41      ROCR_1.0-7           purrr_0.3.4          htmlwidgets_1.5.1   
 [91] patchwork_1.0.0      cowplot_1.0.0        bit_4.0.4            tidyselect_1.0.0     RcppAnnoy_0.0.16    
 [96] magrittr_2.0.1       R6_2.5.0             IRanges_2.16.0       gplots_3.0.4         DBI_1.1.1           
[101] pillar_1.5.1         fitdistrplus_1.1-1   survival_3.2-3       tsne_0.1-3           future.apply_1.5.0  
[106] tibble_3.1.0         crayon_1.4.1         KernSmooth_2.23-18   utf8_1.1.4           plotly_4.9.2.1      
[111] viridis_0.5.1        progress_1.2.2       grid_3.5.1           blob_1.2.1           digest_0.6.27       
[116] tidyr_1.0.2          gridGraphics_0.5-1   stats4_3.5.1         munsell_0.5.0        viridisLite_0.3.0   
[121] ggplotify_0.0.5      

Additional package (directly called in sub-functions):
WGCNA 1.69

```
