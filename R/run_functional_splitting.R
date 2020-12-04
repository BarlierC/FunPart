#' MAIN FUNCTION: functional splitting
#'
#' @param scm single cell matrix (cells in columns / genes in rows)
#' @param tfs vector of TFs names / ID depending on the genes in the scm
#' @param gda BP GO file used for the enrichment analysis
#' @param norm T = perform seurat normalization (default), F = do not perform the normalization
#' @param qtarget top X target genes to consider for each identified cliques 
#' @param adjMethod multiple test correction to perform, either BH (default) or bonferroni
#' @param cutoff p-adjusted value cutoff for the functional enrichment
#' @param percExp percentage of cells a gene needs to be expressed in to be considered
#' @param qExp quantile of gene expression that will be considered as too lowly expressed and so removed
#' @return functional_split object
#' @export
#' @import pheatmap
#' @import data.tree
#' @import plyr
#' @import igraph
#' @import EnvStats
#' @import data.table
#' @import stringr
#' @import clusterProfiler
#' @import Seurat
#' @import Matrix
#' @import doParallel
#' @import foreach
#' @author Celine Barlier
run_functional_splitting <- function(scm,tfs,gda,norm=T,qtarget=0.90,adjMethod="bonferroni",cutoff=0.05,percExp=10,qExp=0.05){ 
    
    ######## INITIALISATION OF THE TREE #########
    
    #Remove ERCC or mt genes if any
    scm <- scm[which(!str_detect(rownames(scm),"^ERCC")),]
    scm <- scm[which(!str_detect(rownames(scm),"^mt")),]
    
    #Not normalized data : necessary for DE GO enrichment
    scm_nm <- scm
    
    #Normalize & clean data
    scm <- filter_norm_data(as.matrix(scm),norm,percExp,qExp)
    
    #Build network & get set of genes based on cliques
    cliqueRes <- find_mutual_cliques(as.matrix(scm),tfs,qtarget)
    
    #If we find cliques
    if(length(cliqueRes) > 0){
      #Perform GO enrichment 
      tmp <- data.frame("SetGenes"=rep(" ",length(cliqueRes)))
      #Structure: gene_set1,gene_set1,gene_set1&gene_set2,gene_set2,gene_set2
      for (i in seq(1,length(cliqueRes))){tmp$SetGenes[i] <- paste(c(paste(unique(c(names(cliqueRes[[i]][[1]]),as.character(c(unlist(cliqueRes[[i]][[1]]))))),collapse=",")),c(paste(unique(c(names(cliqueRes[[i]][[2]]),as.character(c(unlist(cliqueRes[[i]][[2]]))))),collapse=",")),sep="&")}
      
      best_set_genes_go <- runGOanalysis_listSetGenes(as.data.frame(scm_nm),tmp,gda,adjMethod,cutoff)
      
      if(length(best_set_genes_go)>0){
        top_set_genes <- c(strsplit(strsplit(best_set_genes_go$SetGenes[1],split="&")[[1]][1],",")[[1]],strsplit(strsplit(best_set_genes_go$SetGenes[1],split="&")[[1]][2],",")[[1]])
        rm(tmp)
        
        #HC - heatmap
        scm <- as.data.frame(scm)
        mhc <- scm[which(rownames(scm) %in% top_set_genes),]
        
        #If no SD == 0, use it
        if(length(which(colSums(mhc) == 0))==0){
          use <- TRUE
          clNum <- 1
        }else if(length(which(colSums(mhc) == 0))>0 & length(best_set_genes_go$SetGenes) > 1){
          #If SD == 0 but other cliques exists, test them
          rc <- test_new_clique_from_list(scm,best_set_genes_go)
          use <- rc[[1]]
          mhc <- rc[[2]]
          top_set_genes <- rc[[3]]
          clNum <- rc[[4]]
          cl <- rc[[5]]
        }else{
          #If SD == 0 and no other clique - stop
          use <- FALSE
        }
        
        if(use){
          #Clustering + Get first level 
          h <- hierarchical_clust(mhc)
          h <- as.hclust(h)
          cl <- cutree(h,k=2)
          
          #If only one cell is splitted, go to the second level: if more than 5 cells okay, if less try another module. The cells alone is removed = NA
          if(table(cl)[1] == 1 || table(cl)[2] == 1){
            cl <- cutree(h,k=3)
            freqCl <- count(cl)
            na_cell <- names(cl[which(cl == freqCl$x[which(freqCl$freq == min(freqCl$freq))])])
            cl <- cl[which(cl != freqCl$x[which(freqCl$freq == min(freqCl$freq))])]
            if(3 %in% cl){
              cl[cl == 3] <- 2
            }
            
            #Initialization object for output
            iniObj <- init_obj(cl,top_set_genes,best_set_genes_go,cliqueRes)
            
            #Recursive construction of the tree
            t <- build_levels_based_clique(iniObj[[1]],as.matrix(scm),scm_nm,T,c(),2,tfs,iniObj[[2]],iniObj[[3]],gda,iniObj[[4]],adjMethod,cutoff,qtarget,qprobInt,posRatio)
            
            #Object
            functional_res <- build_object(t,scm)
          }else if(table(cl)[1] < 5 || table(cl)[2] < 5){
            #If less than 5 cells in one of the two clusters, try another module if other exists
            rc <- test_new_clique_from_list(scm,best_set_genes_go,clNum,5)
            usethen <- rc[[1]]
            mhc <- rc[[2]]
            top_set_genes <- rc[[3]]
            cl <- rc[[5]]
            if(usethen){
              #Initialization object for output
              iniObj <- init_obj(cl,top_set_genes,best_set_genes_go,cliqueRes)
              
              #Recursive construction of the tree
              t <- build_levels_based_clique(iniObj[[1]],as.matrix(scm),scm_nm,T,c(),2,tfs,iniObj[[2]],iniObj[[3]],gda,iniObj[[4]],adjMethod,cutoff,qtarget,qprobInt,posRatio)
              
              #Object
              functional_res <- build_object(t,scm)
            }else{
              #If nothing is found: empty object
              functional_res <- build_object(list(),scm)
            }
          }else{
            #Initialization object for output
            iniObj <- init_obj(cl,top_set_genes,best_set_genes_go,cliqueRes)
            
            #Recursive construction of the tree
            t <- build_levels_based_clique(iniObj[[1]],as.matrix(scm),scm_nm,T,c(),2,tfs,iniObj[[2]],iniObj[[3]],gda,iniObj[[4]],adjMethod,cutoff,qtarget,qprobInt,posRatio)
            
            #Object
            functional_res <- build_object(t,scm)
          }
        }else{
          #Empty Object
          functional_res <- build_object(list(),scm)
        }
      }else{
        #Empty Object
        functional_res <- build_object(list(),scm)
      }
    }else{
      #Empty Object
      functional_res <- build_object(list(),scm)
    }
    return(functional_res)
  }